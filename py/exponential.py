import numpy as np

mu = 5.0
gamma = 7.0

population_size = 10**6
population = [{'time': 0, 'dead': False} for i in range(population_size)]

for p in population:
    time_to_death = np.random.exponential(scale = mu)
    time_to_growth = np.random.exponential(scale = gamma)

    if time_to_growth > time_to_death:
        p['dead'] = True
        p['time'] = time_to_death
    else:
        p['dead'] = False
        p['time'] = time_to_growth

alive = [p for p in population if p['dead'] == False]
dead = [p for p in population if p['dead'] == True]

alive_percent = 100 * len(alive) / float(population_size)
dead_percent  = 100 * len(dead) / float(population_size)

alive_time = sum([p['time'] for p in population if p['dead'] == False]) / float(len(alive))
dead_time  = sum([p['time'] for p in population if p['dead'] == True]) / float(len(dead))

print '{:10d} alive ({:.2f} percent, average time {:.4f}) \n{:10d} dead  ({:.2f} percent, average time {:.4f}) \n'.format(len(alive), alive_percent, alive_time, len(dead), dead_percent, dead_time)

print((1/gamma) / ((1/mu) + (1/gamma))**2)
