#!/usr/bin/env python3
"""Independent bounded static and fake-interface checks; no R execution or draws."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import socket
import sys
import tempfile
from types import SimpleNamespace
from unittest.mock import patch
import freeze
import run_batch
import summarize

HERE = Path(__file__).resolve().parent

def require_refusal(fn, match):
    try:
        fn()
    except (RuntimeError, ValueError, AssertionError, FileExistsError) as exc:
        assert match in str(exc), str(exc)
    else:
        raise AssertionError('Expected refusal: ' + match)

def static_checks():
    rows, target = freeze.build_campaign_manifest(HERE, HERE)
    recipes = list(csv.DictReader((HERE / 'contract/RECIPES.csv').open()))
    text = (HERE / 'generate.R').read_text()
    specs = dict(re.findall(r'^\s*([\w.]+) = spec\((.*)\),?$', text, re.M))
    assert len(specs) == 64
    inventory = []
    for r in recipes:
        rid = r['recipe_id']; tokens = freeze.split_top_level_commas(specs[rid])
        opts = dict(t.split('=', 1) for t in tokens[4:] if '=' in t)
        opts = {k.strip():v.strip() for k,v in opts.items()}
        assert json.loads(tokens[0]) == r['stream_key']
        assert int(tokens[1].rstrip('L')) == int(r['stream_index'])
        for key, col, default in (('N','n_units',2000), ('T','n_dates',5), ('p','n_covariates',5)):
            assert int(opts.get(key, str(default)).rstrip('L')) == int(r[col]), (rid,key)
        if r['calendar'] != '1:T':
            assert opts['calendar'].replace(' ','') == r['calendar'].replace(' ','')
        for key in ('cohort_law','untreated_law','effect_law','obs_law','error_law','derived'):
            if key in opts:
                value = json.loads(opts[key])
                assert text.count('"' + value + '"') >= 2, (rid,key,value)
        inventory.append({'recipe_id':rid, 'stream_index':int(r['stream_index']), 'explicit_options':opts})
    byid = {r['task_id']: r for r in rows}
    assert len(byid) == len(rows)
    fits = [r for r in rows if r['call_type'] == 'FIT']
    assert len(fits) == 562 and len({r['cell_id'] for r in fits}) == 562
    for row in rows:
        if row['parent_task_id']: assert byid[row['parent_task_id']]['call_type'] == 'FIT'
        if row['call_type'] == 'PRETEST': continue
        assert set(row['primary_coordinate_set'].split(';')) <= set(target['tasks'][row['task_id']]['primary_truth_map'])
    for ct in ('DID_CC', 'DID_UNBALANCED', 'DID_RCS'):
        keys = [(r['recipe_id'], r['literal_call']) for r in rows if r['call_type'] == ct]
        assert len(keys) == len(set(keys)), ct
    # Effective defaults and all important key distinctions are independently asserted.
    args = freeze.default_fit_args('~X1')
    assert freeze.upstream_key(args) == freeze.upstream_key({**args,'weightsname':'NULL','faster_mode':'TRUE'})
    for key, value in (('base_period','"universal"'),('anticipation','1'),('weightsname','"w"'),
                       ('fix_weights','"base_period"'),('faster_mode','FALSE'),('compute_inffunc','FALSE')):
        assert freeze.upstream_key(args) != freeze.upstream_key({**args,key:value}), key
    require_refusal(lambda: freeze.recipe_population_info('UNKNOWN', recipes[0], '', ''), 'Unhandled recipe_id')
    # Independent analytic check for the continuous-cohort linear logit integral.
    p = (math.log1p(math.exp(.8))-math.log1p(math.exp(-.4)))/1.2
    assert abs(freeze.continuous_cohort_mass()-p) < 1e-13 and abs(p-.5) > .01
    # Five prespecified weighted profiles, all with fixed population masses.
    weighted = [r for r in fits if r['primary_target_id'] != 'ORIGINAL_COHORT_ATT']
    assert len(weighted) == 5
    for row in weighted:
        if row['recipe_id'] == 'DATED_WEIGHTS':
            mass = target['tasks'][row['task_id']]['population_masses_map']
            assert mass == {'0':.31, '3':.26, '4':.265, '5':.215}
    dynamic = next(r for r in rows if r['agg_kind'] == 'dynamic')
    meta = target['tasks'][dynamic['task_id']]
    assert all(float(k[4:]) >= 0 for k in meta['primary_truth_map'] if k.startswith('egt:'))
    assert meta['auxiliary_truth_map']
    # The actual R producer exposes the fields consumed by Python, with vector arrays preserved.
    worker = (HERE/'worker.R').read_text()
    for key in ('ci_lower','ci_upper','band_lower','band_upper','overall_ci_lower','overall_ci_upper',
                'ci_lower_egt','ci_upper_egt','band_lower_egt','band_upper_egt','CvM','CvMpval'):
        assert key in worker
    assert 'cells <- lapply(cells, I)' in worker and 'egt = I(as.double(a$egt))' in worker
    assert '1.0.0.0' in worker and '0.1.1.9000' not in worker
    assert 'base::extSoftVersion()' in worker and 'source_identity <- assert_source_identity' in worker
    assert 'md5:' in worker and 'seeds_rds_hash = file_hash' in worker
    assert 'readRDS(task$parent_rds_path)' in worker
    assert 'RNGkind("L\'Ecuyer-CMRG", "Inversion", "Rejection")' in worker
    return inventory

def fake_interfaces(scratch):
    toolkit = scratch/'toolkit'
    shutil.copytree(HERE, toolkit, ignore=shutil.ignore_patterns('__pycache__'))
    # Portability: no external contract/reference/source context is supplied.
    run_batch.subprocess.run([sys.executable, str(toolkit/'freeze.py'), '--check'], check=True, capture_output=True)
    repo = scratch/'checkout'; (repo/'R').mkdir(parents=True); (repo/'.git').mkdir()
    sha = 'a'*40
    (repo/'.git/HEAD').write_text(sha+'\n')
    (repo/'DESCRIPTION').write_text('Package: cdid\nVersion: 1.0.0.0\n')
    (repo/'NAMESPACE').write_text('')
    (repo/'R/example.R').write_text('example <- function(x) x\n')
    assert freeze.checkout_commit(repo) == sha
    fake = scratch/'fake_rscript'
    fake.write_text('''#!/usr/bin/env python3
import sys,json,hashlib
from pathlib import Path
assert sys.argv[1:4][0]=='--vanilla'
assert Path(sys.argv[2]).name=='worker.R'
assert sys.argv[3]=='--freeze'
cfg=json.loads(Path(sys.argv[4]).read_text())
p=Path(cfg['seeds_rds_path']); p.write_bytes(b'FAKE_SEED_ARTIFACT_NO_RNG')
Path(cfg['runtime_out_json']).write_text(json.dumps({'ok':True,'source_sha':cfg['source_sha'],'n_states':29500,
 'seeds_rds_hash':'md5:'+hashlib.md5(p.read_bytes()).hexdigest(),'profile_digest':'FAKE'}))
''')
    fake.chmod(0o700)
    out = scratch/'frozen'
    args = SimpleNamespace(source_sha=sha,out=str(out),repo=str(repo),lib=[str(scratch)],rscript=str(fake),force_new=False)
    rows,targets = freeze.build_campaign_manifest(toolkit,toolkit)
    prov = freeze.run_actual_freeze(args,toolkit,toolkit,rows,targets)
    assert prov['seeds_rds_sha256'] == hashlib.sha256(b'FAKE_SEED_ARTIFACT_NO_RNG').hexdigest()
    run_batch.verify_provenance_on_resume(out,toolkit,toolkit/'manifest.csv',repo_dir=str(repo))
    require_refusal(lambda: freeze.run_actual_freeze(args,toolkit,toolkit,rows,targets), 'empty output')
    seeds=out/'seeds_29500.rds'; saved=seeds.read_bytes();seeds.write_bytes(b'changed')
    require_refusal(lambda: run_batch.verify_provenance_on_resume(out,toolkit,toolkit/'manifest.csv',repo_dir=str(repo)), 'Seed hash')
    seeds.write_bytes(saved)
    ref=toolkit/'reference/CHECKS.json'; saved=ref.read_bytes(); ref.unlink()
    require_refusal(lambda: run_batch.verify_provenance_on_resume(out,toolkit,toolkit/'manifest.csv',repo_dir=str(repo)), 'hash mismatch')
    ref.write_bytes(saved)
    (repo/'R/example.R').write_text('changed')
    require_refusal(lambda: run_batch.verify_provenance_on_resume(out,toolkit,toolkit/'manifest.csv',repo_dir=str(repo)), 'source identity')
    (repo/'R/example.R').write_text('example <- function(x) x\n')
    # Exercise the production launcher argv and real supervisor request/receipt fields.
    launcher = scratch/'fake_task_rscript'
    launcher.write_text("#!/usr/bin/env python3\nimport sys,json\nfrom pathlib import Path\nassert sys.argv[1]=='--vanilla' and Path(sys.argv[2]).name=='worker.R' and sys.argv[3]=='--task'\nt=json.loads(Path(sys.argv[4]).read_text())\nassert t['source_sha']=='a'*40 and Path(t['formals_path']).is_file()\nassert Path(t['repo_path']).joinpath('DESCRIPTION').is_file() and t['frozen_runtime']['profile_digest']=='FAKE'\nassert Path(t['seeds_rds_path']).is_file()\nif t['call_type']=='AGG':\n assert Path(t['parent_rds_path']).read_text()=='FAKE_PARENT'\n result={'agg_result':{'overall_att':1.2,'overall_se':.1,'overall_ci_lower':1.,'overall_ci_upper':1.4,'egt':[],'att_egt':[]}}\nelse:\n assert t['save_parent_rds']\n Path(t['parent_rds_path']).write_text('FAKE_PARENT')\n result={'parent_saved':True,'fit_result':{'cells':{'group':[3],'t':[3],'att':[1.2],'se':[.1],'ci_lower':[1.],'ci_upper':[1.4]}}}\nrec={'attempt_id':t['attempt_id'],'task_id':t['task_id'],'identity_ok':True,'disposition':'success',\n 'startup_seconds':.001,'data_seconds':.001,'operation_seconds':.001,**result}\nPath(t['worker_out_json']).write_text(json.dumps(rec))\n")
    launcher.chmod(0o700)
    receipt = run_batch.run_batch(str(repo), toolkit/'manifest.csv',out,max_calls=2,rscript=str(launcher))
    assert receipt['launched_calls']==2
    results=sorted((out/'records').glob('*.result.json'))
    assert len(results)==2
    records=[json.loads(f.read_text()) for f in results]
    parent=next(r for r in records if r['call_type']=='FIT')
    child=next(r for r in records if r['call_type']=='AGG')
    assert parent['parent_rds_path']==child['parent_rds_path'] and parent['parent_saved']
    assert all(r['disposition']=='success' for r in records)
    assert len(list((out/'batches').glob('*.json')))==1
    # A receipt with an invalid disposition settles as an immutable failure.
    bad=scratch/'bad_worker.py'
    bad.write_text("import sys,json;from pathlib import Path;t=json.loads(Path(sys.argv[-1]).read_text());Path(t['worker_out_json']).write_text(json.dumps({'attempt_id':t['attempt_id'],'task_id':t['task_id'],'identity_ok':True,'disposition':'invalid'}))")
    badout=scratch/'badout'
    first=rows[:1];mf=scratch/'one.csv'
    with mf.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(first[0]));w.writeheader();w.writerows(first)
    with patch.object(run_batch,'expand_manifest_tasks',return_value=run_batch.expand_manifest_tasks(mf,max_reps_override=1)):
        run_batch.run_batch('',mf,badout,allow_unfrozen=True,fake_worker_cmd=[sys.executable,str(bad)])
        failed=list((badout/'records').glob('*.result.json'))
        assert len(failed)==1 and json.loads(failed[0].read_text())['disposition']=='unexpected-error'
        run_batch.run_batch('',mf,badout,allow_unfrozen=True,fake_worker_cmd=[sys.executable,str(bad)])
        assert len(list((badout/'records').glob('*.result.json')))==1
    # No-lock recovery refuses an actually live supervisor and missing/unknown identities.
    orphan=scratch/'orphan';(orphan/'records').mkdir(parents=True)
    opened=orphan/'records/x.open.json'
    opened.write_text(json.dumps({'attempt_id':'x','owner_id':'old','hostname':socket.gethostname(),'supervisor_pid':os.getpid()}))
    require_refusal(lambda: run_batch.acquire_or_reconcile_lock(orphan,'new'), 'alive/unknown')
    # Atomic immutable records survive an attempted overwrite byte-for-byte.
    f=scratch/'immutable.json';run_batch.write_json_atomic({'error':'preserve'},f,overwrite=False); saved=f.read_bytes()
    require_refusal(lambda: run_batch.write_json_atomic({'error':'replace'},f,overwrite=False),'immutable')
    assert f.read_bytes()==saved
    return {'standalone_static_copy':True,'freeze_argv_receipt_hash_resume':True,'missing_seed_source_refusal':True,
            'orphan_liveness_refusal':True,'immutable_failure_record':True,'production_launcher_parent_RDS_paths':True,'malformed_receipt_retained_no_retry':True}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--scratch',required=True)
    args=parser.parse_args(); dest=Path(args.scratch).resolve();dest.mkdir(parents=True,exist_ok=True)
    inventory=static_checks()
    with tempfile.TemporaryDirectory(dir=dest,prefix='interfaces_') as tmp:
        checks=fake_interfaces(Path(tmp))
        # Empty full campaign also exercises every equality pairing/truth key.
        summary=summarize.summarize_campaign(HERE/'manifest.csv', HERE/'targets.json',Path(tmp))
        assert summary['n_attempts']==0
    print(json.dumps({'status':'INDEPENDENT_CHECK_OK','checks':checks,'recipe_inventory':inventory},indent=2))
if __name__=='__main__':main()
