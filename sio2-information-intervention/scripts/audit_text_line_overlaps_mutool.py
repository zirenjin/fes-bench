#!/usr/bin/env python3
"""Supplementary text-line overlap check using MuPDF's structured-text CLI.
This is not a replacement for nature-figure's PyMuPDF collision audit; it checks
only line-box intersections and cannot detect text crossed by graphics.
"""
from pathlib import Path
import argparse, subprocess, tempfile, xml.etree.ElementTree as ET, json

def main():
 ap=argparse.ArgumentParser();ap.add_argument('pdf');ap.add_argument('--out');a=ap.parse_args();pdf=Path(a.pdf)
 with tempfile.TemporaryDirectory(prefix='fes-stext-') as td:
  xml=Path(td)/'text.xml'
  run=subprocess.run(['mutool','draw','-F','stext','-o',str(xml),str(pdf)],capture_output=True,text=True)
  if run.returncode: raise SystemExit(run.stderr or f'mutool exit {run.returncode}')
  root=ET.parse(xml).getroot(); pages=[]; candidates=[]
  for page in root.findall('page'):
   lines=[]
   for block in page.findall('block'):
    for line in block.findall('line'):
     txt=''.join(c.attrib.get('c','') for c in line.findall('.//char')).strip()
     if not txt: continue
     box=[float(v) for v in line.attrib['bbox'].split()]
     lines.append({'text':txt,'bbox':box})
   pages.append({'page':page.attrib.get('id'),'lines':len(lines)})
   for i,x in enumerate(lines):
    for y in lines[i+1:]:
     ix=max(0,min(x['bbox'][2],y['bbox'][2])-max(x['bbox'][0],y['bbox'][0]))
     iy=max(0,min(x['bbox'][3],y['bbox'][3])-max(x['bbox'][1],y['bbox'][1]))
     if ix>=0.5 and iy>=0.5:
      candidates.append({'page':page.attrib.get('id'),'text_a':x['text'],'bbox_a':x['bbox'],'text_b':y['text'],'bbox_b':y['bbox'],'overlap_pt':[ix,iy]})
 result={'pdf':str(pdf),'method':'MuPDF structured text line bbox intersection only','collision_verdict':'NO_TEXT_LINE_OVERLAPS' if not candidates else 'REVIEW_TEXT_LINE_OVERLAPS','threshold_pt':0.5,'pages':pages,'candidate_pairs':candidates,'scope_limit':'Does not inspect text-path intersections, clipping or image/fill overlap; official PyMuPDF audit remains required.'}
 txt=json.dumps(result,indent=2);print(txt)
 if a.out:Path(a.out).write_text(txt+'\n')
 if candidates:raise SystemExit(1)
if __name__=='__main__':main()
