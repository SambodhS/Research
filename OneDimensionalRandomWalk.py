import random
import math
import numpy as np

matrixSize = 1001
iterations = 1000
omega = 0.005
epsilon = 0
k = 0
times = np.logspace(0.1, 4.0, num=20)
for i in range(len(times)):
    times[i] = math.floor(times[i])
print(times)

values1 = []
values2 = []
values3 = []
values4 = []
values5 = []
values6 = []
values7 = []
values8 = []
values9 = []
values10 = []
values11 = []
values12 = []
values13 = []
values14 = []
values15 = []
values16 = []
values17 = []
values18 = []
values19 = []
values20 = []


class App:
    @staticmethod
    def generateMatrix():
        mat = [0] * matrixSize
        mat[math.floor(matrixSize / 2)] = 1
        return mat

    def weights(strength):
        V = omega * strength
        if V > 709:
            V = 709
        weight = np.exp(V)
        return weight

    @staticmethod
    def persistence(persistence, orientation):
        if orientation == "L":
            persistence[1] = (np.exp(-k) / (np.exp(-k) + np.exp(k)))
            persistence[0] = 1 - (np.exp(-k) / (np.exp(-k) + np.exp(k)))
        elif orientation == "R":
            persistence[0] = (np.exp(-k) / (np.exp(-k) + np.exp(k)))
            persistence[1] = 1 - (np.exp(-k) / (np.exp(-k) + np.exp(k)))
        return persistence

    def memoryReduction(mem):
        i = 0
        while i < len(mem):
            mem[i] = mem[i] - epsilon * mem[i]
            if mem[i] < 0:
                mem[i] = 0
            i += 1
        return mem

    @staticmethod
    def getOrientation():
        r = random.uniform(0, 1)
        if r <= 0.5:
            orientation = "L"
        else:
            orientation = "R"
        return orientation

    def moveAgent(mat, agentPosition, pL, pR, wL, wR, persistence, orientation):
        r = random.uniform(0, 1)
        if r < pL:
            if agentPosition != 0:
                mat[agentPosition - 1] = 1
                agentPosition = agentPosition - 1
                orientation = "L"
            else:
                print("Oh no!")
                App.moveAgent(mat, agentPosition, 0, 1, wL, wR, persistence, orientation)
        elif pL < r <= pL + pR:
            if agentPosition != matrixSize - 1:
                mat[agentPosition + 1] = 1
                agentPosition = agentPosition + 1
                orientation = "R"
            else:
                print("Oh no!")
                App.moveAgent(mat, agentPosition, 1, 0, wL, wR, persistence, orientation)
        return mat, agentPosition, orientation

    @staticmethod
    def main():
        for iteration in range(iterations):
            print(iteration)
            mat = App.generateMatrix()
            mem = App.generateMatrix()
            persistence = [0] * 2
            orientation = App.getOrientation()
            agentPosition = int(matrixSize / 2)
            for steps in range(times[19]):
                try:
                    wL = App.weights(mem[agentPosition - 1])
                except IndexError:
                    wL = 0
                try:
                    wR = App.weights(mem[agentPosition + 1])
                except IndexError:
                    wR = 0
                persistence = App.persistence(persistence, orientation)
                mat[agentPosition] = 0
                mat, agentPosition, orientation = App.moveAgent(mat, agentPosition, ((wL * persistence[0]) / ((wL * persistence[0]) + (wR * persistence[1]))), ((wR * persistence[1]) / ((wL * persistence[0]) + (wR * persistence[1]))), wL, wR, persistence, orientation)
                mem[agentPosition] = mem[agentPosition] + 1
                #mem = App.memoryReduction(mem)
                if steps == times[0] - 1:
                    values1.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[1] - 1:
                    values2.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[2] - 1:
                    values3.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[3] - 1:
                    values4.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[4] - 1:
                    values5.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[5] - 1:
                    values6.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[6] - 1:
                    values7.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[7] - 1:
                    values8.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[8] - 1:
                    values9.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[9] - 1:
                    values10.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[10] - 1:
                    values11.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[11] - 1:
                    values12.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[12] - 1:
                    values13.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[13] - 1:
                    values14.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[14] - 1:
                    values15.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[15] - 1:
                    values16.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[16] - 1:
                    values17.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[17] - 1:
                    values18.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[18] - 1:
                    values19.append(agentPosition - math.floor(matrixSize / 2))
                elif steps == times[19] - 1:
                    values20.append(agentPosition - math.floor(matrixSize / 2))
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values1[w] ** 2)
            mstd = mstd + (values1[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[0], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values2[w] ** 2)
            mstd = mstd + (values2[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[1], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values3[w] ** 2)
            mstd = mstd + (values3[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[2], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values4[w] ** 2)
            mstd = mstd + (values4[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[3], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values5[w] ** 2)
            mstd = mstd + (values5[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[4], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values6[w] ** 2)
            mstd = mstd + (values6[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[5], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values7[w] ** 2)
            mstd = mstd + (values7[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[6], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values8[w] ** 2)
            mstd = mstd + (values8[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[7], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values9[w] ** 2)
            mstd = mstd + (values9[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[8], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values10[w] ** 2)
            mstd = mstd + (values10[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[9], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values11[w] ** 2)
            mstd = mstd + (values11[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[10], msd, merr)
        msd = 0.0
        mstd = 0.0
        for w in range(iterations):
            msd = msd + (values12[w] ** 2)
            mstd = mstd + (values12[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[11], msd, merr)
        for w in range(iterations):
            msd = msd + (values13[w] ** 2)
            mstd = mstd + (values13[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[12], msd, merr)
        for w in range(iterations):
            msd = msd + (values14[w] ** 2)
            mstd = mstd + (values14[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[13], msd, merr)
        for w in range(iterations):
            msd = msd + (values15[w] ** 2)
            mstd = mstd + (values15[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[14], msd, merr)
        for w in range(iterations):
            msd = msd + (values16[w] ** 2)
            mstd = mstd + (values16[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[15], msd, merr)
        for w in range(iterations):
            msd = msd + (values17[w] ** 2)
            mstd = mstd + (values17[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[16], msd, merr)
        for w in range(iterations):
            msd = msd + (values18[w] ** 2)
            mstd = mstd + (values18[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[17], msd, merr)
        for w in range(iterations):
            msd = msd + (values19[w] ** 2)
            mstd = mstd + (values19[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[18], msd, merr)
        for w in range(iterations):
            msd = msd + (values20[w] ** 2)
            mstd = mstd + (values20[w] ** 2) ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[19], msd, merr)


if __name__ == "__main__":
    App.main()
