import numpy as np

mu = 0.9
gamma = 0.05
sigma = 0.05

population_size = 10**6
population = [{'stage': 0, 'age': 0, 'dead': False} for i in range(population_size)]

while 0 in [p['stage'] for p in population]:
    for p in population:
        if p['stage'] == 0:
            fate = np.random.multinomial(1, [mu, gamma, sigma])
            if fate[0]:
                p['dead'] = True
                p['stage'] = -1
                p['age'] += 1
            elif fate[1]:
                p['stage'] += 1
                p['age'] += 1
            elif fate[2]:
                p['age'] += 1

alive = [p for p in population if p['dead'] == False]
dead = [p for p in population if p['dead'] == True]

alive_percent = 100 * len(alive) / float(population_size)
dead_percent  = 100 * len(dead) / float(population_size)

alive_time = sum([p['age'] for p in population if p['stage'] == 1]) / float(len(alive))
dead_time  = sum([p['age'] for p in population if p['stage'] == -1]) / float(len(dead))

print '{:10d} alive ({:.2f} percent, average age {:.4f}) \n{:10d} dead  ({:.2f} percent, average age {:.4f}) \n'.format(len(alive), alive_percent, alive_time, len(dead), dead_percent, dead_time)

print(1.0 / (mu + gamma))
