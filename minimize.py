from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize,shgo,differential_evolution

E = 210e9  # E-Modul in Pa
rho = 7850  # Dichte in kg/m^
price = 10

#use cross-entropy method to optimize.
def cem(num_steps,iterations,batch_size,waning_time= 10,additional_std=0.25,fraction=0.2):
    #mean rewards
    mean_rewards = []
    # initialize mean and standard deviation
    theta_mean = np.array([7.5,7.5,1.25,1.25,0.0325])
    cons_cov = np.array([2.5,2.5,0.75,0.75,0.03])*0.1
    theta_cov = np.array([2.5,2.5,0.75,0.75,0.03])*0.1

    for iteration in range(iterations):
        #initialize  batch
        theta_set=[]
        rewards=[]
        for batch in range(batch_size):
            #for each batch
            #X = AZ + u
            #theta=(theta_cov + (np.identity(5)*max(1 - iteration / waning_time, 0) * additional_std**2)) @ np.random.randn(5) +theta_mean
            theta = (theta_cov + (cons_cov * max(1 - iteration / waning_time, 0) * additional_std ** 2)) @ np.random.randn(5) + theta_mean
            print(theta)
            theta_set.append(theta)
            reward = cost(theta)
            if tension(theta) >=0.2e9:
                continue
            rewards.append(reward)
        # process that selects elite sets
        index = np.argsort(np.array(rewards))
        index = index[::-1]
        elite = []
        num_top = int(fraction * batch_size)
        for i in range (num_top):
            elite.append(theta_set[index[i]])
        matrix =np.stack(elite,axis=0).T
        for i in range(5):
            theta_mean[i] = np.mean(matrix[i,:])

        # cov = AA.T, using eigenvalue decomposition
        e,v = np.linalg.eig(np.cov(matrix))
        for i in range(5):
            if e[i] <0 : e[i] =0

        e = np.identity(5) * np.sqrt(e)
        theta_cov = v @ e

        mean_rewards.append(np.mean(rewards))
    return mean_rewards,theta_mean

def tension(x: np.ndarray):
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    for i in range(-1, -5, -1):
        myCrane.addExternalForce(i, 0, 0, -500e3 / 4)
    fem = FEM(myCrane)
    # N, R, U = fem.TrussAnalysis()
    # t = np.max(N[np.newaxis]) / x[4]
    t = np.max(fem.getTension())
    print(t * 1e-9)
    return -t


def cost(x: np.ndarray):
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    return np.sum(myCrane.mass)


if __name__ == "__main__":
    #print(tension(np.array([7.5, 7.5, 1, 1, 0.0225])))
    cons = (
        {'type': 'ineq', 'fun': lambda x: tension(x) + 0.2e9},  # tension < 0.2e9
    )
    res = shgo(cost,
                   bounds=((5, 10), (5, 10), (0.5, 2), (0.5, 2), (2.5e-3, 6.25e-2)),
                   constraints=cons)
    print(res)
    #rewards,x=cem(100,30,30)
    #plt.plot(rewards)
    #plt.legend()
    #plt.show()
    #[ 7.500e+00  7.500e+00  1.250e+00  1.250e+00  3.250e-02]