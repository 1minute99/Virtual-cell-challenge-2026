"""Real systemd integration checks: isolation, literal argv, duplicate guard, exit status."""
import json
import subprocess
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LAUNCHER = ROOT / 'scripts/run_detached.sh'

class DetachedRunTest(unittest.TestCase):
    def test_detachment_and_duplicate_protection(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'result.json'
            name = f'recovery-test-{time.time_ns()}'
            unit = f'aivc-{name}.service'
            payload = 'literal $(touch BAD) `echo BAD` $HOME'
            program = 'import os,sys,time,json,pathlib; time.sleep(1); pathlib.Path(sys.argv[1]).write_text(json.dumps([os.getppid(),os.isatty(1),sys.argv[2]]))'
            args = ['bash',str(LAUNCHER),name,directory,'/usr/bin/python3','-c',program,str(out),payload]
            try:
                first = subprocess.run(args,capture_output=True,text=True)
                self.assertEqual(first.returncode,0, first.stderr)
                duplicate = subprocess.run(args,capture_output=True,text=True)
                self.assertNotEqual(duplicate.returncode,0)
                for _ in range(100):
                    if out.exists():break
                    time.sleep(.1)
                self.assertTrue(out.exists(),first.stdout)
                parent, tty, actual = json.loads(out.read_text())
                self.assertFalse(tty)
                self.assertEqual(actual,payload)
                self.assertEqual(Path(f'/proc/{parent}/comm').read_text().strip(), 'systemd')
                result = subprocess.check_output(['systemctl','--user','show',unit,'-p','Result','-p','Slice','-p','ExecMainStatus'],text=True)
                self.assertIn('Result=success',result)
                self.assertIn('ExecMainStatus=0',result)
                self.assertIn('Slice=aivc-jobs.slice',result)
            finally:
                subprocess.run(['systemctl','--user','stop',unit],capture_output=True)
                subprocess.run(['systemctl','--user','reset-failed',unit],capture_output=True)

    def test_failure_is_visible(self):
        with tempfile.TemporaryDirectory() as directory:
            name=f'recovery-failure-{time.time_ns()}'
            unit=f'aivc-{name}.service'
            try:
                launch=subprocess.run(['bash',str(LAUNCHER),name,directory,'/bin/sh','-c','exit 7'],capture_output=True,text=True)
                self.assertEqual(launch.returncode,0,launch.stderr)
                for _ in range(100):
                    result=subprocess.check_output(['systemctl','--user','show',unit,'-p','Result','-p','ExecMainStatus'],text=True)
                    if 'ExecMainStatus=7' in result:break
                    time.sleep(.1)
                self.assertIn('ExecMainStatus=7',result)
                self.assertIn('Result=exit-code',result)
            finally:
                subprocess.run(['systemctl','--user','stop',unit],capture_output=True)
                subprocess.run(['systemctl','--user','reset-failed',unit],capture_output=True)

if __name__ == '__main__':unittest.main()
