const {chromium}=require('[local archive path]');
const fs=require('fs'),path=require('path');
(async()=>{
 const root=process.argv[2]||path.resolve(__dirname,'..');
 const browser=await chromium.launch({executablePath:'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',headless:true,args:['--disable-gpu']});
 const page=await browser.newPage({viewport:{width:1440,height:1100}}),errors=[],missing=[];
 page.on('pageerror',e=>errors.push(e.message));page.on('requestfailed',r=>missing.push(r.url()));
 await page.goto('file:///'+root.replace(/\\/g,'/')+'/index.html#probe-labels',{waitUntil:'load'});
 try { await page.waitForFunction(()=>document.documentElement.dataset.probeReady==='true'&&document.documentElement.dataset.ready==='english'); } catch(e) { console.error(JSON.stringify({errors,missing,globals:await page.evaluate(()=>({probe:!!window.FES_PROBE_LABELS,data:!!window.FES_DATA,ready:document.documentElement.dataset.ready,probeReady:document.documentElement.dataset.probeReady,body:document.body.innerText.slice(0,250)}))}));throw e; }
 await page.evaluate(()=>document.documentElement.style.scrollBehavior='auto');
 const result={counts:await page.evaluate(()=>window.FES_PROBE_LABELS.counts)};
 const count=()=>page.locator('#probe-labels tbody tr:visible').count();
 if(await count()!==25)throw Error('Expected 25 phase-layer rows');
 result.systems={};for(const s of ['SiO2','Ti','Zr','Hf','CaSiO3','MgSiO3','Fe','C']){await page.selectOption('#probe-system',s);result.systems[s]=await count();if(!result.systems[s])throw Error('Empty system '+s);}
 await page.selectOption('#probe-system','all');await page.selectOption('#probe-view','thermal');result.scopedThermal=await count();if(result.scopedThermal!==5)throw Error('Expected Ti plus 4 carbon scoped thermal rows');
 await page.selectOption('#probe-view','blocked');result.blocked=await count();if(result.blocked!==3)throw Error('Expected Ca cubic, Mg PPv and historical graphite');
 await page.selectOption('#probe-view','all');await page.fill('#probe-search','C2221');result.c2221=await count();if(result.c2221!==1)throw Error('C2221 duplicate or missing');
 await page.fill('#probe-search','no-such-reference-xyz');if(await count()!==0||!await page.locator('#probe-label-detail').innerText().then(t=>t.includes('No rows')))throw Error('Empty filter broken');
 await page.fill('#probe-search','');await page.locator('[data-probe-select="PL14"]').first().click();
 const detail=await page.locator('#probe-label-detail').innerText();if(!detail.includes('30 GPa'))throw Error('Expected Mg Pv pressure gate');
 const sources=page.locator('#probe-label-detail .probe-source button');await sources.first().click();if(!await page.locator('#source-dialog').isVisible())throw Error('Source evidence does not open');await page.locator('#source-close').click();
 for(const id of await page.locator('#probe-labels tbody tr').evaluateAll(es=>es.map(e=>e.dataset.probeId))){await page.locator('[data-probe-select="'+id+'"]').first().click();if(await page.locator('#probe-label-detail').getAttribute('data-row-id')!==id)throw Error('Detail mismatch '+id);}
 await page.locator('[data-probe-select="PL01"]').first().click();await page.locator('#probe-labels').evaluate(e=>e.scrollIntoView({block:'start',behavior:'instant'}));
 await page.screenshot({path:path.join(root,'probe-matrix-desktop-preview.png')});
 result.downloads=await page.locator('#probe-labels .downloads a').evaluateAll(es=>es.map(e=>e.getAttribute('href')));
 for(const href of result.downloads){if(!fs.existsSync(path.join(root,href)))throw Error('Missing download '+href);}
 await page.setViewportSize({width:390,height:844});await page.locator('#probe-labels').evaluate(e=>e.scrollIntoView({block:'start',behavior:'instant'}));await page.screenshot({path:path.join(root,'probe-matrix-mobile-preview.png')});
 result.mobile=await page.evaluate(()=>({width:innerWidth,body:document.body.scrollWidth,table:document.querySelector('.probe-table').getBoundingClientRect().width}));
 if(result.mobile.body>result.mobile.width+1)throw Error('Mobile overflow');
 result.errors=errors;result.missing=missing;fs.writeFileSync(path.join(root,'PROBE_BROWSER_CHECK.json'),JSON.stringify(result,null,2));await browser.close();console.log(JSON.stringify(result,null,2));if(errors.length||missing.length)process.exitCode=1;
})().catch(e=>{console.error(e);process.exit(1)});
