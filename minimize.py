from crane import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize,shgo,differential_evolution

E = 210e9  # E-Modul in Pa
rho = 7850  # Dichte in kg/m^
price = 10

#use cross-entropy method to optimize.
def cem(iterations,batch_size,waning_time= 20,additional_std=0.45,fraction=0.4):
    #mean rewards
    mean_rewards = []
    # initialize mean and standard deviation
    theta_mean = np.array([7.5,7.5,1.25,1.25,0.0325])
    cons_cov = np.diagonal([[2.5,2.5,0.75,0.75,0.03]])
    A= np.identity(5)

    for iteration in range(iterations):
        #initialize  batch
        theta_set=[]
        rewards=[]
        for batch in range(batch_size):
            #for each batch
            #X = AZ + u
            #theta=(theta_cov + (np.identity(5)*max(1 - iteration / waning_time, 0) * additional_std**2)) @ np.random.randn(5) +theta_mean
            theta = (A + (cons_cov * max(1 - iteration / waning_time, 0) * additional_std ** 2)) @ np.random.randn(5) + theta_mean
            theta=np.abs(theta)
            theta[0] = np.clip(theta[0], 5, 10)
            theta[1] = np.clip(theta[1], 5, 10)
            theta[2] = np.clip(theta[2],0.5,2)
            theta[3] = np.clip(theta[3],0.5,2)
            theta[4] = np.clip(theta[4],0.025,0.0625)
            if tension(theta) >=0.2e9:
                continue
            theta_set.append(theta)
            reward = cost(theta)
            rewards.append(reward)
        # process that selects elite sets
        index = np.argsort(np.array(rewards))
        #reverse sequence, so we have [987654321]
        index = index[::-1]
        elite = []
        num_top = int(fraction * len(theta_set))
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
        A = v @ e
        mean_rewards.append(np.mean(rewards))
    return mean_rewards,theta_mean

def tension(x: np.ndarray):
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    for i in range(-1, -5, -1):
        myCrane.addExternalForce(i, 0, 0, -500e3 / 4)
    fem = FEM(myCrane)
    # N, R, U = fem.TrussAnalysis()
    # t = np.max(N[np.newaxis]) / x[4]
    t = np.abs(np.max(fem.getTension()))
    return t


def cost(x: np.ndarray):
    myCrane = crane(1, x[0], x[1], x[2], x[3], x[4], rho, E)
    return -np.sum(myCrane.mass)


if __name__ == "__main__":
    #print(tension(np.array([7.5, 7.5, 1, 1, 0.0225])))
    # cons = (
    #     {'type': 'ineq', 'fun': lambda x: tension(x) + 0.2e9},  # tension < 0.2e9
    # )
    # res = shgo(cost,
    #                bounds=((5, 10), (5, 10), (0.5, 2), (0.5, 2), (2.5e-3, 6.25e-2)),
    #                constraints=cons)
    # print(res)
    rewards,x=cem(100,30)
    plt.plot(rewards,x)
    plt.legend()
    plt.show()