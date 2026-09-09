"""Calibre ru_maxrss contre VmHWM après lancement d'un interpréteur frais.

Diagnostic Linux : ne mesure aucun solveur et ne classe pas leur mémoire.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys


def controler():
    code = ('import json,resource; '
            'print(json.dumps(dict(rusage_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'
            'proc=[s.strip() for s in open("/proc/self/status") '
            'if s.startswith(("VmHWM:","VmRSS:"))])))')
    observations = []
    for octets in (0, 160*1024*1024):
        reserve = bytearray(octets)
        r = subprocess.run([sys.executable, "-c", code], stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, check=True)
        observations.append(dict(parent_reserve_octets=len(reserve),
            parent_rusage_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            enfant=json.loads(r.stdout)))
    return dict(python=sys.version, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                observations=observations,
                portee="Diagnostic du compteur, hors campagne ; aucune mesure de mémoire des solveurs")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sortie", type=Path)
    a = p.parse_args()
    a.sortie.write_text(json.dumps(controler(), ensure_ascii=False, indent=2)+"\n")
