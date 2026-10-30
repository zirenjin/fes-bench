"""Build a provenance-bound, read-only DFT explorer from existing archives.

No VASP/MD/model calculations are submitted. Gamma eigenvectors are derived
from archived force constants, and checked against archived Gamma frequencies.
Licensed POTCAR bodies, wavefunctions, and charge densities are not exported.
"""
from pathlib import Path
import csv, json, hashlib, re, shutil, collections, os, xml.etree.ElementTree as ET
import numpy as np
import yaml

HERE=Path(__file__).resolve().parent
PROJECT=Path(os.environ.get('FES_ARCHIVE_PROJECT_ROOT', '.')).expanduser().resolve()
BASE=PROJECT/'outputs/fes-bench-reference-v2'
RESET=PROJECT/'outputs/fes-label-gap-reset-20261021'
DATA=HERE/'data'; DATA.mkdir(exist_ok=True)
PROVENANCE=[]; ERRORS=[]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def source(p):
    p=Path(p); digest=sha(p)
    try: path=str(p.relative_to(PROJECT))
    except ValueError: path=p.name
    item={'archive_path':path,'sha256':digest,'bytes':p.stat().st_size}
    if item not in PROVENANCE: PROVENANCE.append(item)
    return item

def public_value(x):
    if isinstance(x,dict):return {k:public_value(v) for k,v in x.items()}
    if isinstance(x,list):return [public_value(v) for v in x]
    if isinstance(x,str):
        x=re.sub(r'/mnt/[a-z]/[^\s\";]*?/dpa-adapt/', 'project/', x, flags=re.I)
        x=re.sub(r'/expanse/lustre/projects/[^/]+/[^/]+/', 'Expanse/project/', x)
        x=re.sub(r'/home/[^/]+/', 'Expanse/home/', x)
        x=re.sub(r'/GenSIvePFS/users/[^/]+/', 'GenSi/project/', x)
    return x

def dump(path,obj): path.write_text(json.dumps(public_value(obj),separators=(',',':'),allow_nan=False)+'\n')
def emit(name,obj):
    dump(DATA/(name+'.json'),obj)
    (DATA/(name+'.js')).write_text('window.DFT_DATA=window.DFT_DATA||{};window.DFT_DATA['+json.dumps(name)+']='+json.dumps(public_value(obj),separators=(',',':'),allow_nan=False)+';\n')
def rows(p): return list(csv.DictReader(open(p)))
def numeric(x):
    try:
        v=float(x); return v if np.isfinite(v) else None
    except (TypeError,ValueError): return None
def arr(e):
    return [[float(x) for x in (v.text or '').split()] for v in e.findall('v')] if e is not None else []
def clean(a,decimals=7): return np.round(np.asarray(a),decimals).tolist()
def infer(p):
    s=str(p).lower()
    system=next((v for k,v in [('casio3','CaSiO3'),('ca-tetragonal','CaSiO3'),('ca_nvt','CaSiO3'),('mgsio3','MgSiO3'),('sio2','SiO2'),('carbon','C'),('c_dft','C'),('ti-hcp','Ti'),('hf','Hf'),('/ti/','Ti'),('/zr/','Zr'),('/fe/','Fe'),('w6-fe','Fe')] if k in s),'Other')
    phase=next((v for k,v in [('graphite','graphite'),('diamond','diamond'),('post_perovskite','PPv'),('/ppv/','PPv'),('/pv/','Pv'),('tetragonal','tetragonal'),('cubic','cubic'),('quartz','quartz'),('cristobalite','cristobalite'),('c2221','C2221'),('tridymite','tridymite'),('hcp','hcp'),('bcc','bcc'),('fcc','fcc')] if k in s),'unspecified')
    if 'ca_nvt' in s:phase='cubic candidate'
    return system,phase

def xml_record(p):
    root=ET.parse(p).getroot()
    species=[(r.find('c').text or '').strip() for r in root.findall('./atominfo/array[@name="atoms"]/set/rc')]
    pars={i.get('name'):(i.text or '').strip() for i in root.findall('./incar/*') if i.get('name')}
    system,phase=infer(p)
    if system=='Other':
        elements=set(species)
        for sy,els in [('C',{'C'}),('Ti',{'Ti'}),('Zr',{'Zr'}),('Hf',{'Hf'}),('Fe',{'Fe'}),('MgSiO3',{'Mg','Si','O'}),('CaSiO3',{'Ca','Si','O'}),('SiO2',{'Si','O'})]:
            if elements==els:system=sy;break
    frames=[];scf=[]
    for k,c in enumerate(root.findall('calculation')):
        structure=c.find('structure'); positions=arr(structure.find('varray[@name="positions"]')) if structure is not None else []
        cell=arr(structure.find('crystal/varray[@name="basis"]')) if structure is not None else []
        force=arr(c.find('varray[@name="forces"]'))
        energy={i.get('name'):numeric(i.text) for i in c.findall('energy/i')}
        se=[{i.get('name'):numeric(i.text) for i in step.findall('energy/i')} for step in c.findall('scstep')]
        scf.extend({'ionic_step':k+1,'electronic_step':j+1,'energy_eV':x.get('e_fr_energy')} for j,x in enumerate(se) if x.get('e_fr_energy') is not None)
        if positions and cell:
            frames.append({'step':k+1,'cell':cell,'fractional':positions,'forces':force,'energy':energy,'electronic_steps':len(se),'volume_A3':float(abs(np.linalg.det(cell))),'max_force_eV_A':float(np.linalg.norm(force,axis=1).max()) if force else None})
    if not frames:raise ValueError('No complete ionic structure in XML')
    footer=False;out=p.with_name('OUTCAR')
    if out.exists():
        with out.open('rb') as f:f.seek(max(0,out.stat().st_size-16000));footer=b'General timing and accounting' in f.read()
    rec={'id':'vasp-'+sha(p)[:16],'system':system,'phase':phase,'name':str(p.parent.relative_to(PROJECT/'outputs')),'kind':'DFT ionic/force calculation','n_atoms':len(species),'symbols':species,'parameters':pars,'frames':frames,'scf':scf,'normal_footer':footer,'source':source(p),'frame_count':len(frames),'status':'normal footer present' if footer else 'footer not verified','admission':'Raw output availability is separate from label admission.'}
    return rec

def xdat_record(p):
    lines=p.read_text(errors='replace').splitlines();scale=float(lines[1].split()[0]);cell=np.array([[float(x) for x in s.split()[:3]] for s in lines[2:5]])*scale
    sy=lines[5].split();counts=[int(x) for x in lines[6].split()];n=sum(counts);symbols=sum(([s]*c for s,c in zip(sy,counts)),[]);frames=[]
    i=7
    while i<len(lines):
        if 'configuration' in lines[i].lower():
            try:pos=np.array([[float(x) for x in s.split()[:3]] for s in lines[i+1:i+1+n]])
            except ValueError:break
            if pos.shape==(n,3):frames.append({'step':len(frames)+1,'cell':clean(cell),'fractional':clean(pos),'forces':[],'energy':{},'volume_A3':float(abs(np.linalg.det(cell))),'max_force_eV_A':None,'electronic_steps':None})
            i+=n+1
        else:i+=1
    if not frames:raise ValueError('No XDATCAR frames')
    system,phase=infer(p);meta={'id':'traj-'+sha(p)[:16],'system':system,'phase':phase,'name':str(p.parent.relative_to(PROJECT/'outputs')),'kind':'DFT MD trajectory','n_atoms':n,'symbols':symbols,'parameters':{},'frames':frames,'scf':[],'normal_footer':None,'source':source(p),'frame_count':len(frames),'status':'archived coordinates','admission':'Recorded dynamics; not an accepted free-energy label. Frame index is not a convergence certificate.'}
    out=p.with_name('OUTCAR')
    if out.exists():
        text=out.read_text(errors='replace')
        times=re.findall(r'POTIM\s*=\s*([\d.E+\-]+)',text[:120000]);meta['POTIM_fs']=float(times[0]) if times else None
        energies=[float(x) for x in re.findall(r'FREE ENERGIE OF THE ION-ELECTRON SYSTEM.*?free\s+energy\s+TOTEN\s*=\s*([\d.E+\-]+)',text,re.S)]
        blocks=re.findall(r'POSITION\s+TOTAL-FORCE \(eV/Angst\)[^\n]*\n\s*-{10,}[^\n]*\n(.*?)(?=^\s*-{10,})',text,re.S|re.M)
        forces=[];out_positions=[]
        for block in blocks:
            vals=[line.split() for line in block.strip().splitlines()]
            try:
                forces.append([[float(x) for x in row[3:6]] for row in vals if len(row)==6])
                out_positions.append([[float(x) for x in row[:3]] for row in vals if len(row)==6])
            except ValueError:forces.append([]);out_positions.append([])
        # Bind arrays only when a one-to-one ionic correspondence is available.
        if len(energies)==len(frames):
            for fr,e in zip(frames,energies):fr['energy']={'e_fr_energy':e}
        if len(forces)==len(frames) and all(len(x)==n for x in forces):
            residuals=[]
            for fr,positions in zip(frames,out_positions):
                delta=np.asarray(positions)@np.linalg.inv(cell)-np.asarray(fr['fractional']);delta-=np.round(delta);residuals.append(float(np.linalg.norm(delta@cell,axis=1).max()))
            meta['max_position_alignment_error_A']=max(residuals)
            if max(residuals)<1e-3:
                for fr,f in zip(frames,forces):fr['forces']=f;fr['max_force_eV_A']=float(np.linalg.norm(f,axis=1).max())
        meta['normal_footer']='General timing and accounting' in text[-16000:];meta['outcar_source']=source(out)
        meta['force_frame_alignment']='one-to-one bound' if all(len(fr['forces'])==n for fr in frames) else 'unavailable; coordinates only'
    return meta

def fc_array(p):
    with p.open() as f:
        shape=[int(x) for x in f.readline().split()];n=shape[0];m=shape[-1];fc=np.zeros((n,m,3,3));row_ids=[]
        for block in range(n*m):
            line=f.readline()
            while line and not line.strip():line=f.readline()
            ij=[int(x)-1 for x in line.split()]
            if block%m==0:row_ids.append(ij[0])
            fc[block//m,ij[1]]=np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
    return fc,np.array(row_ids)

def phonon_record(p,geometry=None,fcfile=None):
    y=geometry or yaml.load(p.read_text(),Loader=yaml.CSafeLoader if hasattr(yaml,'CSafeLoader') else yaml.SafeLoader);pr=y.get('primitive_cell') or y.get('unit_cell');su=y.get('supercell')
    if not pr or not su:raise ValueError('Missing primitive/supercell geometry')
    fp=fcfile or p.with_name('FORCE_CONSTANTS')
    if not fp.exists():raise ValueError('Missing matching FORCE_CONSTANTS')
    fc,row_ids=fc_array(fp);lat=np.array(pr['lattice']);frac=np.array([x['coordinates'] for x in pr['points']]);mass=np.array([x['mass'] for x in pr['points']]);spos=np.array([x['coordinates'] for x in su['points']])@np.array(su['lattice'])@np.linalg.inv(lat)
    mapping=[]
    for point in spos:
        delta=point-frac;delta-=np.round(delta);dist=np.linalg.norm(delta@lat,axis=1);j=int(dist.argmin())
        if dist[j]>1e-4:raise ValueError('Supercell/primitive mapping mismatch')
        mapping.append(j)
    mapping=np.array(mapping);n=len(frac);D=np.zeros((n,3,n,3))
    for a in range(n):
        i=int(np.flatnonzero(mapping[row_ids]==a)[0])
        for b in range(n):D[a,:,b,:]=fc[i,mapping==b].sum(axis=0)/np.sqrt(mass[a]*mass[b])
    D=D.reshape(3*n,3*n);anti=float(np.linalg.norm(D-D.T));vals,vec=np.linalg.eigh((D+D.T)/2);factor=float(y.get('phonopy',{}).get('frequency_unit_conversion_factor',15.633302));freq=np.sign(vals)*np.sqrt(abs(vals))*factor
    qfile=p.with_name('qpoints.yaml');qerr=None
    if not qfile.exists():qfile=p.with_name('mesh.yaml')
    if qfile.exists():
        q=yaml.load(qfile.read_text(),Loader=yaml.CSafeLoader if hasattr(yaml,'CSafeLoader') else yaml.SafeLoader);q0=next((a for a in q.get('phonon',[]) if np.linalg.norm(a.get('q-position',[1,1,1]))<1e-8),None)
        if q0 and len(q0['band'])==len(freq):qerr=float(np.max(abs(np.sort(freq)-np.sort([a['frequency'] for a in q0['band']]))))
    disp=vec.T.reshape(3*n,n,3)/np.sqrt(mass)[None,:,None]
    disp/=np.linalg.norm(disp,axis=2).max(axis=1)[:,None,None]
    dataset={'id':'phonon-'+sha(fp)[:16],'system':infer(p)[0],'phase':infer(p)[1],'name':str(p.parent.relative_to(PROJECT/'outputs')),'cell':clean(lat),'fractional':clean(frac),'symbols':[x['symbol'] for x in pr['points']],'frequencies_THz':clean(freq),'displacements':clean(disp),'gamma_frequency_check_max_abs_THz':qerr,'dynamical_asymmetry_norm':anti,'source':source(fp),'geometry_source':source(p),'supercell_atoms':len(mapping),'ifc_shape':list(fc.shape[:2]),'ifc_block_norm_eV_A2':clean(np.linalg.norm(fc,axis=(2,3)),5),'response_first_atom_per_A':clean(-fc[0].transpose(1,0,2)),'response_convention':'Uses Phi(i,0)=Phi(0,i)^T reciprocity for compact force constants.','supercell_cell':su['lattice'],'supercell_fractional':[x['coordinates'] for x in su['points']],'supercell_symbols':[x['symbol'] for x in su['points']],'mode_scope':'Gamma-point eigenvectors derived from archived force constants. Oscillation amplitudes are visual, not thermal. Negative modes are displacement illustrations, not stable oscillators.'}
    if qerr is not None and qerr>1e-3:raise ValueError(f'Gamma frequency comparison failed: {qerr}')
    mesh=p.with_name('mesh.yaml')
    if mesh.exists():
        # C-backed YAML loader avoids slow pure-Python parsing of full meshes.
        my=yaml.load(mesh.read_text(),Loader=yaml.CSafeLoader if hasattr(yaml,'CSafeLoader') else yaml.SafeLoader)
        frequencies=[];weights=[]
        for q in my.get('phonon',[]):
            for band in q['band']:frequencies.append(band['frequency']);weights.append(q['weight'])
        edges=np.linspace(min(-1.,min(frequencies)),max(frequencies)+.01,121);hist,edges=np.histogram(frequencies,bins=edges,weights=weights,density=False)
        dataset['dos']={'frequency_THz':clean((edges[1:]+edges[:-1])/2),'mode_weight':hist.tolist(),'bin_width_THz':float(edges[1]-edges[0]),'q_weight_sum':sum(q['weight'] for q in my['phonon']),'negative_mode_weight':float(np.sum(np.array(weights)[np.array(frequencies)<-1e-3])),'source':source(mesh),'normalization':'Raw q-multiplicity-weighted mode counts per frequency bin; not divided by atom count.'}
    dataset['status']='Historical superseded FC reconstruction; not a thermodynamic reference' if 'phonon_assembly_v05_retry3' in str(p) else 'Archived force-constant dataset; thermodynamic admission is separate'
    return dataset

def ti_geometry(p,super_path):
    def parse(path):
        a=path.read_text().splitlines();scale=float(a[1]);lat=np.array([[float(x) for x in l.split()[:3]] for l in a[2:5]])*scale
        try:counts=[int(x) for x in a[5].split()];k=6
        except ValueError:counts=[int(x) for x in a[6].split()];k=7
        if a[k].lower().startswith('s'):k+=1
        direct=a[k].lower().startswith('d');k+=1;pos=np.array([[float(x) for x in l.split()[:3]] for l in a[k:k+sum(counts)]])
        if not direct:pos=pos*scale@np.linalg.inv(lat)
        return {'lattice':lat.tolist(),'points':[{'symbol':'Ti','coordinates':v.tolist(),'mass':47.867} for v in pos]}
    source(super_path)
    return {'primitive_cell':parse(p),'supercell':parse(super_path),'phonopy':{'frequency_unit_conversion_factor':15.633302}}

def main():
    ledger=rows(BASE/'COMPUTE_LEDGER.csv');ledger_src=source(BASE/'COMPUTE_LEDGER.csv')
    refresh_path=DATA/'remote_refresh.json';refresh=json.load(open(refresh_path)) if refresh_path.exists() else None
    if refresh:
        ledger+=refresh['new_ledger_records']
    ledger=public_value(ledger)
    for x in ledger:
        # Keep archive state verbatim, including failures and non-DFT work.
        s=' '.join([x.get('task_id',''),x.get('vasp',''),x.get('failure_or_notes','')]).lower()
        x['work_class']='postprocess / reference / model sampling' if any(k in x['task_id'].lower() for k in ['einstein','spring-scan','postprocess','nvt-smoke','nvt-continuation','dp-gauss','dp-midpoint']) and 'vasp' not in x['task_id'].lower() else 'DFT or DFT workflow record'
    inventory=[];seen={};xmls=sorted(set(BASE.rglob('vasprun.xml'))|set(RESET.rglob('vasprun.xml')))
    for p in xmls:
        digest=sha(p)
        if digest in seen:
            seen[digest]['archive_copies'].append(str(p.relative_to(PROJECT)));continue
        try:
            r=xml_record(p);r['archive_copies']=[str(p.relative_to(PROJECT))];emit(r['id'],r)
            brief={k:r[k] for k in ['id','system','phase','name','kind','n_atoms','parameters','normal_footer','source','frame_count','status','archive_copies']};inventory.append(brief);seen[digest]=brief
        except Exception as ex:ERRORS.append({'path':str(p.relative_to(PROJECT)),'error':str(ex),'kind':'XML parse'})
    remote_path=DATA/'remote_raw_export.json';remote=json.load(open(remote_path)) if remote_path.exists() else None
    remote_added=0
    if remote:
        for r in remote['records']:
            digest=r['source']['sha256']
            if digest in seen:continue
            emit(r['id'],r);brief={k:r[k] for k in ['id','system','phase','name','kind','n_atoms','parameters','normal_footer','source','frame_count','status','archive_copies']};inventory.append(brief);seen[digest]=brief;PROVENANCE.append(r['source']);remote_added+=1
        ERRORS.extend(remote['errors'])
    # Full merged Ca trajectory plus separate bracket/extension records. These
    # are linked records, never summed as independent runs or samples.
    ca=BASE/'ca_nvt_bracket_analysis_20261024'
    for sub in ['V319_merged_20261027_verified/merged_work','V319_order_ess_extension_20261026','V323_live_0757UTC']:
        p=ca/sub/'XDATCAR'
        if p.exists():
            try:
                r=xdat_record(p);emit(r['id'],r);inventory.insert(0,{k:r[k] for k in ['id','system','phase','name','kind','n_atoms','parameters','normal_footer','source','frame_count','status']})
            except Exception as ex:ERRORS.append({'path':str(p.relative_to(PROJECT)),'error':str(ex),'kind':'trajectory parse'})
    phonons=[];fcseen=set()
    geometry_files=sorted(set(BASE.rglob('phonopy.yaml'))|set(BASE.rglob('phonopy_params.yaml')))
    for p in geometry_files:
        if not p.with_name('FORCE_CONSTANTS').exists():continue
        digest=sha(p.with_name('FORCE_CONSTANTS'))
        if digest in fcseen:continue
        try:
            r=phonon_record(p);emit(r['id'],r);phonons.append({k:r[k] for k in ['id','system','phase','name','supercell_atoms','source','gamma_frequency_check_max_abs_THz','status']});fcseen.add(digest)
        except Exception as ex:ERRORS.append({'path':str(p.relative_to(PROJECT)),'error':str(ex),'kind':'phonon assembly'})
    tiroot=BASE/'w9-ti-hcp-qha-volume-pilot-20261021/remote_archive'
    for fp in sorted(tiroot.glob('*/FORCE_CONSTANTS')):
        p=fp.with_name('POSCAR.unitcell');sp=fp.parent/'perfect/POSCAR'
        if not p.exists() or not sp.exists() or sha(fp) in fcseen:continue
        try:
            r=phonon_record(p,geometry=ti_geometry(p,sp));emit(r['id'],r);phonons.append({k:r[k] for k in ['id','system','phase','name','supercell_atoms','source','gamma_frequency_check_max_abs_THz','status']});fcseen.add(sha(fp))
        except Exception as ex:ERRORS.append({'path':str(p.relative_to(PROJECT)),'error':str(ex),'kind':'Ti phonon assembly'})
    cf=BASE/'w7-carbon-source-audit-20261020'
    fvt=rows(cf/'carbon_p5_qha_candidate_20261025/FVT_CANDIDATE.csv');ext=rows(cf/'carbon_extended_volumes_20261025/EXTENDED_FVT_CANDIDATE.csv')
    therm=[]
    for layer,rs in [('Carbon v1 · 5 GPa bracket',fvt),('Carbon v2 · extended graphite bracket',ext)]:
        for r in rs:
            therm.append({'layer':layer,'system':'C','phase':r['phase'],'T_K':float(r['T_K']),'V_A3_atom':float(r['V_A3_atom']),'E_static_eV_atom':float(r['E_sigma_to_zero_eV_atom']),'F_vib_eV_atom':float(r['F_vib_eV_atom']),'F0_eV_atom':float(r['F0_eV_atom']),'volume_ratio':float(r['V_over_V0'])})
    for p in [cf/'carbon_p5_qha_candidate_20261025/FVT_CANDIDATE.csv',cf/'carbon_extended_volumes_20261025/EXTENDED_FVT_CANDIDATE.csv']:source(p)
    tp=BASE/'w9-ti-hcp-qha-volume-pilot-20261021/remote_qmesh_convergence/TI_HCP_FVT_QHA_CANDIDATE.csv'
    # kJ/mol of Ti2 -> eV per atom, not eV per cell of 36 atoms.
    for r in rows(tp):
        if r['mesh']!='(15, 15, 15)':continue
        therm.append({'layer':'Ti hcp · mesh 15 · harmonic comparison','system':'Ti','phase':'hcp','T_K':float(r['T_K']),'V_A3_atom':float(r['V_A3_Ti2'])/2,'E_static_eV_atom':float(r['E0_kJmol_Ti2'])/96.48533212331002/2,'F_vib_eV_atom':float(r['Fvib_kJmol_Ti2'])/96.48533212331002/2,'F0_eV_atom':float(r['F0_kJmol_Ti2'])/96.48533212331002/2,'volume_ratio':float(r['V_over_V0'])})
    source(tp)
    master=rows(BASE/'MASTER_FREE_ENERGY_GRID_MANIFEST.csv');source(BASE/'MASTER_FREE_ENERGY_GRID_MANIFEST.csv')
    labelrecords=public_value([r for r in master if not r['record_type'].startswith('source_') and not r['record_type'].startswith('derived_source_')])
    phasegrid=[]
    for r in master:
        if r['system']=='C' and r['record_type']=='required_paired_DeltaG' and numeric(r['label_value']) is not None:
            phasegrid.append({'T_K':numeric(r['T_K']),'P_GPa':numeric(r['P_GPa']),'deltaG_graphite_minus_diamond_meV_atom':float(r['label_value'])*1000,'admission':r['admission_state'],'reference_layer':r['reference_layer'],'numerical_error':r['convergence_error'],'notes':r['notes']})
    gates=RESET/'execution_update_20261021/SIX_PHASE_PHONON_ACCEPTANCE_20261022.csv'
    six=rows(gates) if gates.exists() else [];source(gates) if gates.exists() else None
    entropy_path=PROJECT/'outputs/fes-bench-sio2-information-intervention-20261027/data/D2_CARBON_HARMONIC_PROBE_LABELS.csv'
    entropy=rows(entropy_path) if entropy_path.exists() else [];source(entropy_path) if entropy else None
    summary={'date':'2026-10-30','scope':'Archived project calculations plus a timestamped read-only Expanse refresh when present. This is not a continuously live queue. External published labels are not counted as new DFT runs.','ledger_source':ledger_src,'ledger_records':len(ledger),'local_xml_files_scanned':len(xmls),'raw_xml_files_scanned':len(xmls)+(len(remote['records'])+len(remote['errors']) if remote else 0),'unique_xml_records':len(seen),'remote_xml_added':remote_added,'ionic_playback_records':len(inventory),'phonon_records':len(phonons),'ledger_system_counts':dict(collections.Counter(x['system'] or 'Unspecified' for x in ledger)),'raw_system_counts':dict(collections.Counter(x['system'] for x in inventory)),'label_record_counts':dict(collections.Counter(x['system'] for x in labelrecords)),'remote_refresh':refresh,'errors':ERRORS,'capabilities':{'charge_density':'No charge-density export supplied. No electronic isosurface is fabricated.','band_structure':'No dedicated band path is exported. Gamma mode stems are not a phonon dispersion.','live_queue':'Timestamped read-only refresh only; not continuously live.'}}
    if refresh:summary['remote_refresh']={k:v for k,v in refresh.items() if k!='new_ledger_records'}
    emit('catalog',{'summary':summary,'calculations':inventory,'ledger':ledger,'phonons':phonons,'thermal':therm,'carbon_grid':phasegrid,'labels':labelrecords,'six_phase_acceptance':six,'entropy':entropy})
    dump(DATA/'provenance.json',PROVENANCE);dump(DATA/'extraction_errors.json',ERRORS)
    with (DATA/'COMPUTE_LEDGER.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
    with (DATA/'FREE_ENERGY_ADMISSION.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(labelrecords[0]));w.writeheader();w.writerows(labelrecords)
    (DATA/'BUILD_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['ledger_source','errors','capabilities']},indent=2));print('Extraction issues',len(ERRORS),ERRORS[:3])

if __name__=='__main__':main()
