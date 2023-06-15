import crane
from fem import FEM
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
import os

E = 210e9  # E-module [Pa]
rho = 7850  # density [kg/m^3]
price = 10

def calculateCov(range_min:np.ndarray,range_max: np.ndarray):
    """
    This function returns the covariance matrix
    @param range_min:
    @param range_max:
    @return: covariance matrix
    """
    # Calculate variances for each dimension
    variances = ((range_max - range_min)) / 6
    cov = np.diag(variances)
    return cov

def eigDecomposition(cov):
    """
    This function decomposes eigenvector
    @param cov:
    @return:Matrix with eigenvector
    """
    e, v = np.linalg.eig(cov)
    for i in range(4):
        if e[i] < 0: e[i] = 0
    e = np.identity(4) * np.sqrt(e)
    A = v @ e
    return A


def cem(crane,iterations,batch_size,waning_time= 3,additional_std=0.35,fraction=0.15):
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



if __name__ == "__main__":
    result_1=[]
    result_2=[]
    result_3=[]
    t1=[]
    t2=[]
    t3=[]
    cpun=4
    N = cpun*1
    # we use multiprocessors here
    p = Pool(os.cpu_count())
    #resule_1 look like this : [[rewards=[],mean=[]],...,]
    for i in range(N):
        p.apply_async(cem,   args=(crane.crane_1,100,500),callback=result_1.append)
        # p.apply_async(cem, args=(crane.crane_2_1,100,500),callback=result_2.append)
        # p.apply_async(cem, args=(crane.crane_2_2,100,500),callback=result_3.append)
    p.close()
    p.join()
    # for i in range(50):
    #     rewards_1, x_1 = cem(crane.crane_1,100, 500)
    #     rewards_2, x_2 = cem(crane.crane_2_1, 100, 500)
    #     rewards_3, x_3 = cem(crane.crane_2_2, 100, 500)
    #     result_1.append(rewards_1[-1])
    #     result_2.append(rewards_2[-1])
    #     result_3.append(rewards_3[-1])
    #     t1.append(x_1[-1])
    #     t2.append(x_2[-1])
    #     t3.append(x_3[-1])
    # print(t1)
    # print(t2)
    # print(t3)
    rewards_1 = [result_1[i][0] for i in range(N)]
    # rewards_2 = [result_2[i][0] for i in range(N)]
    # rewards_3 = [result_3[i][0] for i in range(N)]
    t1 = np.array([result_1[i][1] for i in range(N)])
    # t2 = np.array([result_2[i][1] for i in range(N)])
    # t3 = np.array([result_3[i][1] for i in range(N)])
    fig, axs = plt.subplots(6, 1, figsize=(8, 6))
    axs[0].plot(rewards_1)
    axs[0].set_title('crane_1_rewards')

    axs[1].scatter(t1[:,-2], t1[:,-1])
    axs[1].set_title('crane_1_volumn')

    # axs[2].plot(rewards_2)
    # axs[2].set_title('crane_2_1_rewards')
    #
    # axs[3].scatter(t2[:,-2], t2[:,-1])
    # axs[3].set_title('crane_2_1_volumn')
    #
    # axs[4].plot(rewards_3)
    # axs[4].set_title('crane_2_2_rewards')
    #
    # axs[5].scatter(t3[:,-2], t3[:,-1])
    # axs[5].set_title('crane_2_2_volumn')

    # rewards,x=cem(100,450)
    # vol=[]
    # for s in x:
    #     vol.append(s[2]*s[3])
    # indexNorm= np.argmin(norm)
    # print("theta: "+ str(x[-1]))
    # print("reward: " + str(rewards[-1]))
    # fig, axs = plt.subplots(2, 1, figsize=(8, 6))
    # axs[0].plot(rewards[1:], label='r')
    # axs[0].set_title('r')
    # axs[0].legend()
    # axs[1].plot(vol, label='volumn')
    # axs[1].set_title('x')
    # axs[1].legend()
    plt.show()

