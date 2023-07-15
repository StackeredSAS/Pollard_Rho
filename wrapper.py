#!/usr/bin/env python3
import subprocess
from secrets import randbelow
import time
import sys
from multiprocessing import cpu_count

# read the order from the header file
# this is the simplest method I found so the user doesn't need to modify 2 files or use a config file
with open("pollard_rho_worker.h", "r") as f:
    for line in f.readlines():
        if line.startswith("#define ORDER"):
            order = int(line.split('"')[1].split('"')[0], 16)

DPOINTS = {}

# choose how many parallel processes you want to launch (default: number of CPU)
num_processes = cpu_count()

def newP():
    return subprocess.Popen(["./worker", str(randbelow(2**64))], stdout=subprocess.PIPE)

def parseStdout(p):
    # U.x,a,b
    stdout, _ = process.communicate()
    lines = stdout.decode().strip().split("\n")
    if len(lines) == 1:
        return []
    return [[int(e, 16) for e in x.split(",")] for x in lines]

def checkFound(p):
    for line in parseStdout(process):
        x, a, b = line
        known = DPOINTS.get(x)
        if known is not None:
            a_, b_ = known
            if b_ != b:
                priv = (a-a_)*pow(b_-b, -1, order)
                priv %= order
                print(f"Found private key ! x = {priv}")
                return True
        else:
            DPOINTS[x] = (a, b)

if __name__ == '__main__':
    processes = [newP() for _ in range(num_processes)]

    while True:
        for i, process in enumerate(processes):
            if process.poll() == 0:
                # Respawn a new process because the previous one has finished
                processes[i] = newP()
                if checkFound(process):
                    # kill all previous processes and exit
                    for remaining_process in processes:
                        if remaining_process.poll() is None:
                            remaining_process.terminate()
                    sys.exit(0)

        # sleep to not consume CPU power endlessly polling
        time.sleep(1)

