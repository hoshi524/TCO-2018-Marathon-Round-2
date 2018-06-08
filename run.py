import subprocess
import json
import threading
import queue
import json
import sys
import os

subprocess.call('g++ --std=c++0x -W -Wall -O2 -s -pipe -mmmx -msse -msse2 -msse3 -o out/main.out {}'.format(sys.argv[1]), shell=True)
subprocess.call('javac -d out src/CrystalLightingVis.java', shell=True)

scorefile = "best-score.json"
case = 20
scores = [0]

try:
    with open(scorefile, "r") as f:
        scores = json.loads(f.read())
except FileNotFoundError:
    scores = [-10000] * case

if case > len(scores):
    scores = [-10000] * case

def solve(seed):
    return float(subprocess.check_output('java -cp out CrystalLightingVis -exec out/main.out -seed {}'.format(seed), shell=True))

class State:
    count = 0
    rate = 0
    lock = threading.Lock()

    def add(self, seed, score):
        if scores[seed] < score:
            scores[seed] = score
        nom = 0
        if scores[seed] > 0:
            nom = score / scores[seed] 
        with self.lock:
            self.count = self.count + 1
            self.rate = (self.rate * (self.count - 1) + nom) / self.count
            print('{}\t{:8.2f}\t{}\t{}'.format(seed, score, nom, self.rate))

state = State()
q = queue.Queue()

def worker():
    while True:
        seed = q.get()
        if seed is None:
            break
        try:
            score = solve(seed)
            state.add(seed, score)
            q.task_done()
        except ValueError as err:
            print("seed {} : {}".format(seed, err))
            os._exit(1)

num_worker_threads = 4
threads = []
for i in range(num_worker_threads):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)

for seed in range(1, case):
    q.put(seed)

q.join()

for i in range(num_worker_threads):
    q.put(None)
for t in threads:
    t.join()

with open(scorefile, "w") as f:
    f.write(json.dumps(scores))