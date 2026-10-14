"""Copy only the stopped output tree to a fresh, authorized data-drive folder."""
from pathlib import Path
from datetime import datetime,timezone
import json,os,subprocess,time

ROOT=Path('/home/sji/Documents/Programming/Python/embryogenesis-simulation')
SOURCE=ROOT/'output';TARGET=Path('/mnt/data/Documents/Python/embryogenesis-simulation/output')
HERE=Path(__file__).parent
if __name__=='__main__':
    stopped=json.loads((HERE/'interruption.json').read_text())
    assert stopped['worker_exited'] and not Path(f'/proc/{stopped["old_pid"]}').exists()
    assert SOURCE.is_dir() and not SOURCE.is_symlink() and not TARGET.exists()
    assert not (HERE/'initial-copy.json').exists()
    TARGET.mkdir()
    argv=['rsync','-aH','--stats',str(SOURCE)+'/',str(TARGET)+'/'];start=time.monotonic()
    with (HERE/'initial-copy.log').open('x') as stream:
        result=subprocess.run(argv,stdout=stream,stderr=subprocess.STDOUT)
    proof=dict(source=str(SOURCE),destination=str(TARGET),created_destination_empty=True,
        initial_rsync_exit_code=result.returncode,seconds=time.monotonic()-start,argv=argv,
        utc=datetime.now(timezone.utc).isoformat())
    (HERE/'initial-copy.json').write_text(json.dumps(proof,indent=2)+'\n')
    assert result.returncode==0,'Stable output tree copy failed; original retained'
    print(json.dumps(proof,indent=2))
