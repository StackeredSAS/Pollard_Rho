#!/usr/bin/env python3
import subprocess
import random
import time
import sys
from multiprocessing import cpu_count

# choose how many parallel processes you want to launch (default: number of CPU)
num_processes = cpu_count()

def newP():
    return subprocess.Popen(["./pollard_rho", str(random.randint(0, 2**64))], stdout=subprocess.PIPE)

if __name__ == '__main__':
    processes = [newP() for _ in range(num_processes)]

    while True:
        for i, process in enumerate(processes):
            return_code = process.poll()

            if return_code == 0:
                stdout, _ = process.communicate()
                print(stdout.decode().strip())
                for remaining_process in processes:
                    if remaining_process.poll() is None:
                        remaining_process.terminate()
                sys.exit(0)

            elif return_code is not None:
                # Respawn a new process if the previous one has finished
                processes[i] = newP()
        # sleep to not consume CPU power endlessly polling
        time.sleep(10)

