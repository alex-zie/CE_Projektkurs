import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import os

E = 210e9  # E-module [Pa]
rho = 7850  # density [kg/m^3]
price = 1

def calculateCov(range_min:np.ndarray,range_max: np.ndarray):
    """
    This function returns the covariance matrix
    @param range_min:
    @param range_max:
    @return: covariance matrix
    """
    # Calculate variances for each dimension
    variances = ((range_max - range_min)) /4
    cov = np.diag(variances)
    return cov

def eigDecomposition(cov):
    """
    This function decomposes eigenvector
    @param cov:
    @return:Matrix with eigenvector
    """
    e, v = np.linalg.eig(cov)
    for i in range(2):
        if e[i] < 0: e[i] = 0
    e = np.identity(2) * np.sqrt(e)
    A = v @ e
    return A


def cem(crane,iterations,batch_size,waning_time= 1,additional_std=0.35,fraction=0.15):
    """
    This function runs cross-entropy method to search the best combination of parameter set that has minimal cost
    @param iterations: total iterations that the function runs
    @param batch_size: total batches in one iteration
    @param waning_time: first <waning_time> iteraitons that additional standard deviation added
    @param additional_std: additional standard deviation
    @param fraction: percentage of the top tier
    @return:
    """
    #set mean rewards
    mean_rewards = []
    mean_rewards.append(0)
    # initialize mean and standard deviation
    theta_mean = np.array([10,10,1.25,0.0325])
    #length and height are set around 10
    A = eigDecomposition(calculateCov(np.array([9.9, 9.9, 0.5, 0.0025]), np.array([10.1, 10.1, 2, 0.0625])))
    cons_cov = eigDecomposition(calculateCov(np.array([9.9, 9.9, 0.5, 0.0025]), np.array([10.1, 10.1, 2, 0.0625])))
    mean_x=[]
    for iteration in range(iterations):
        # initialize batch
        theta_set=[]
        rewards=[]
        print("-----------------------iteration: ",iteration,"--------------------------")
        np.random.seed(None)
        for batch in range(batch_size):
            # for each batch
            # X = AZ + u
            # see multi-variant distribution
            # theta=(theta_cov + (np.identity(5)*max(1 - iteration / waning_time, 0) * additional_std**2)) @ np.random.randn(5) +theta_mean
            x = (A + (cons_cov * max(1 - iteration / waning_time, 0) * additional_std ** 2)) @ np.random.randn(4) + theta_mean
            # normal distribution generate some samples with negative value, need absolute value
            x=np.abs(x)
            # initialize crane object
            # manually selects the samples within certain range
            x[0] = np.clip(x[0], 9.9, 10)
            x[1] = np.clip(x[1], 9.9, 10)
            x[2] = np.clip(x[2],0.5,2)
            x[3] = np.clip(x[3],0.0025,0.0625)
            myCrane = crane( x[0], x[1], x[2], x[3],  rho, E, False)
            for i in range(-1, -5, -1):
                myCrane.addExternalForce(i, 0, 0, -500e3 / 4)
            fem = FEM(myCrane, False)
            # calculate tension and cost
            t = np.abs(np.max(fem.getTension()))
            cost = -np.sum(myCrane.mass)*price
            # check condition
            if t >=0.2e9 or not fem.check_bending_force():
                continue
            # store sample
            theta_set.append(x)
            rewards.append(cost)
        # process that selects elite sets
        index = np.argsort(np.array(rewards))
        #reverse sequence, so we have [987654321]
        index = index[::-1]
        elite = []
        # number of candidates
        num_top = int(fraction * len(theta_set))
        # append candidates
        for i in range (num_top):
            elite.append(theta_set[index[i]])
        # append corresponding parameter sets of candidates
        matrix =np.stack(elite,axis=0).T
        # calculate mean value of all parameter in this elite sets
        for i in range(4):
            theta_mean[i] = np.mean(matrix[i,:])
        # deep copy into returned lists
        mean_x.append(np.copy(theta_mean))
        # generate new covariance for the next iteration
        # cov = AA.T, using eigenvalue decomposition
        A = eigDecomposition(np.cov(matrix))
        diff = np.abs(mean_rewards[-1]-np.mean(rewards))
        print("---------------diff: ",diff,"---------------------")
        if(diff<10):
            break
        mean_rewards.append(np.mean(rewards))
    return [mean_rewards[-1] ,mean_x[-1]]


def sample(crane, mean, A, cons_cov, iteration, waning_time, additional_std):
    np.random.seed(None)
    x = (A + (cons_cov * max(1 - iteration / waning_time, 0) * additional_std ** 2)) @ np.random.randn(2) + mean
    x[1] = np.clip(x[1], 0.5, 1.414)
    x[0] = np.clip(np.abs(x[0]), 0, 10*500e3)
    myCrane = crane(10, 10, x[1], 0.0025, rho, E, False)
    for i in myCrane.tip_nodes:
        myCrane.addExternalForce(i, 0, 0, -500e3/len(myCrane.tip_nodes))
    # counterweight
    for i in myCrane.counterweight_nodes:
        myCrane.addExternalForce(i, 0, 0, -x[0] / len(myCrane.counterweight_nodes))
    fem = FEM(myCrane, own_weight=True)
    fem.optimize_crossections(625e-4, 0.9*200e6)
    cost = -np.sum(fem.truss.mass)
    print("sample: ",x,"after： ",cost)
    return [x, cost]

def cem_slcw(crane,iterations,batch_size,waning_time= 1,additional_std=0.3,fraction=0.15):
    """
    This function runs cross-entropy method to search the best combination of parameter set that has minimal cost
    @param crane: crane version
    @param iterations: total iterations that the function runs
    @param batch_size: total batches in one iteration
    @param waning_time: first <waning_time> iteraitons that additional standard deviation added
    @param additional_std: additional standard deviation
    @param fraction: percentage of the top tier
    @return:
    """
    #set mean rewards
    mean_rewards = []
    mean_rewards.append(0)
    # initialize mean and standard deviation
    theta_mean = np.array([2500e3,0.957])
    #length and height are set around 10
    A = calculateCov(np.array([0, 0.5]), np.array([10*500e3,1.414]))
    mean_x=[]
    for iteration in range(iterations):
        # initialize batch
        results=[]
        print("-----------------------iteration: ",iteration,"--------------------------")
        def transfer(result):
            results.append(result)
            return
        """
        Pool(enter your number of physical cpu)
        """
        p = Pool(os.cpu_count())
        for batch in range(batch_size):
            """
            try to parallelize the loop
            """
            p.apply_async(sample, args=(crane,theta_mean,A,A,iteration,waning_time,additional_std),callback=transfer)

        # process that selects elite sets
        p.close()
        p.join()
        rewards = [results[i][1] for i in range(len(results))]
        theta_set = [results[i][0] for i in range(len(results))]
        index = np.argsort(np.array(rewards))
        #reverse sequence, so we have [987654321]
        index = index[::-1]
        elite = []
        # number of candidates
        num_top = int(fraction * len(theta_set))
        # append candidates
        for i in range (num_top):
            elite.append(theta_set[index[i]])
        # append corresponding parameter sets of candidates
        matrix =np.stack(elite,axis=0).T
        # calculate mean value of all parameter in this elite sets
        for i in range(2):
            theta_mean[i] = np.mean(matrix[i,:])
        # deep copy into returned lists
        mean_x.append(np.copy(theta_mean))
        # generate new covariance for the next iteration
        # cov = AA.T, using eigenvalue decomposition
        A = eigDecomposition(np.cov(matrix))
        diff = np.abs(mean_rewards[-1]-np.mean(rewards))
        print("---------------diff: ",diff,"---------------------")
        if(diff<5):
            break
        mean_rewards.append(np.mean(rewards))
    return mean_rewards ,mean_x




if __name__ == "__main__":
    """
    Print result for one optimization
    """
    rewards, x = cem_slcw(crane.crane_2_1, 100, 1000,waning_time=3,additional_std=1,fraction=0.25)
    print("theta: "+ str(x[-1]))
    print("reward: " + str(rewards[-1]))

    """
    Print result in chronological order
    """
    plt.figure(constrained_layout=True)
    plt.plot(rewards[1:])
    plt.xlabel("Iteration")
    plt.ylabel("Cost")
    plt.show()
    """
    multiple times optimization
    """
    # fig, axs = plt.subplots(2, 1, figsize=(8, 6))
    # cost=[]
    # xset=[]
    # for i in range(10):
    #     rewards, x = cem_slcw(crane.crane_2_2, 100, 1000)
    #     cost.append(rewards[-1])
    #     xset.append(x[-1])
    # axs[0].plot(cost, label='r')
    # axs[0].set_title('rewards')
    # axs[0].legend()
    # axs[1].scatter([xset[i][0] for i in range(len(xset))], [xset[i][1] for i in range(len(xset))])
    # axs[1].set_title('cw and sl')
    # axs[1].legend()
    # plt.show()
