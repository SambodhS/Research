import random
import math
import numpy as np

matrixSize = 501
iterations = 2000
omega = 0.05
epsilon = 0
k = 0

times = np.logspace(0.1, 4.0, num=20)
for i in range(len(times)):
    times[i] = math.floor(times[i])

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
        mat = [[0] * matrixSize for _ in range(matrixSize)]
        mat[math.floor(matrixSize / 2)][math.floor(matrixSize / 2)] = 1
        return mat

    def getAgentPosition(mat):
        agentPosition = [0] * 2
        i = 0
        while i < len(mat):
            j = 0
            while j < len(mat[i]):
                if mat[i][j] == 1:
                    agentPosition[0] = i
                    agentPosition[1] = j
                j += 1
            i += 1
        return agentPosition

    def weights(strength):
        V = omega * strength
        if V > 709:
            V = 709
        weight = np.exp(V)
        return weight

    @staticmethod
    def getOrientation():
        r = random.uniform(0, 1)
        if r <= 0.25:
            orientation = "U"
        elif 0.25 < r <= 0.5:
            orientation = "D"
        elif 0.5 < r <= 0.75:
            orientation = "L"
        else:
            orientation = "R"
        return orientation

    @staticmethod
    def persistence(persistence, orientation):
        if orientation == "U":
            persistence[0] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif orientation == "D":
            persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[1] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif orientation == "L":
            persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[2] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif orientation == "R":
            persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            persistence[3] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        return persistence

    def memoryReduction(mem):
        i = 0
        while i < len(mem):
            j = 0
            while j < len(mem[i]):
                if mem[i][j] > 0:
                    mem[i][j] = mem[i][j] - epsilon * mem[i][j]
                j += 1
            i += 1
        return mem

    def moveAgent(mat, agentPosition, pU, pD, pL, pR, wU, wD, wL, wR, persistence, orientation):
        r = random.uniform(0, 1)
        if r < pU:
            if agentPosition[0] != 0:
                mat[agentPosition[0] - 1][agentPosition[1]] = 1
                agentPosition[0] = agentPosition[0] - 1
                orientation = "U"
            else:
                print("Oh no!", agentPosition)
                App.moveAgent(mat, agentPosition, 0, (wD * persistence[1]) / (wD * persistence[1] + wL * persistence[2] + wR * persistence[3]), (wL * persistence[2]) / (wD * persistence[1] + wL * persistence[2] + wR * persistence[3]), (wR * persistence[3]) / (wD * persistence[1] + wL * persistence[2] + wR * persistence[3]), wU, wD, wL, wR, persistence, orientation)
        elif pU < r <= pU + pD:
            if agentPosition[0] != matrixSize - 1:
                mat[agentPosition[0] + 1][agentPosition[1]] = 1
                agentPosition[0] = agentPosition[0] + 1
                orientation = "D"
            else:
                print("Oh no!", agentPosition)
                App.moveAgent(mat, agentPosition, (wU * persistence[0]) / (wU * persistence[0] + wL * persistence[2] + wR * persistence[3]), 0, (wL * persistence[2]) / (wU * persistence[0] + wL * persistence[2] + wR * persistence[3]), (wR * persistence[3]) / (wU * persistence[0] + wL * persistence[2] + wR * persistence[3]), wU, wD, wL, wR, persistence, orientation)
        elif pU + pD < r <= pU + pD + pL:
            if agentPosition[1] != 0:
                mat[agentPosition[0]][agentPosition[1] - 1] = 1
                agentPosition[1] = agentPosition[1] - 1
                orientation = "L"
            else:
                print("Oh no!", agentPosition)
                App.moveAgent(mat, agentPosition, (wU * persistence[0]) / (wU * persistence[0] + wD * persistence[1] + wR * persistence[3]), (wD * persistence[1]) / (wU * persistence[0] + wD * persistence[1] + wR * persistence[3]), 0, (wR * persistence[3]) / (wU * persistence[0] + wD * persistence[1] + wR * persistence[3]), wU, wD, wL, wR, persistence, orientation)
        elif pU + pD + pL < r <= pU + pD + pL + pR:
            if agentPosition[1] != matrixSize - 1:
                mat[agentPosition[0]][agentPosition[1] + 1] = 1
                agentPosition[1] = agentPosition[1] + 1
                orientation = "R"
            else:
                print("Oh no!", agentPosition)
                App.moveAgent(mat, agentPosition, (wU * persistence[0]) / (wU * persistence[0] + wD * persistence[1] + wL * persistence[2]), (wD * persistence[1]) / (wU * persistence[0] + wD * persistence[1] + wL * persistence[2]), (wL * persistence[2]) / (wU * persistence[0] + wD * persistence[1] + wL * persistence[2]), 0, wU, wD, wL, wR, persistence, orientation)
        return mat, agentPosition, orientation

    def print2D(mat):
        for row in mat:
            print(row)

    @staticmethod
    def main():
        for iteration in range(iterations):
            print(iteration)
            mat = App.generateMatrix()
            mem = App.generateMatrix()
            persistence = [0] * 4
            orientation = App.getOrientation()
            agentPosition = [math.floor(matrixSize / 2), math.floor(matrixSize / 2)]
            for steps in range(int(times[19])):
                try:
                    wU = App.weights(mem[agentPosition[0] - 1][agentPosition[1]])
                except IndexError:
                    wU = 0
                try:
                    wD = App.weights(mem[agentPosition[0] + 1][agentPosition[1]])
                except IndexError:
                    wD = 0
                try:
                    wL = App.weights(mem[agentPosition[0]][agentPosition[1] - 1])
                except IndexError:
                    wL = 0
                try:
                    wR = App.weights(mem[agentPosition[0]][agentPosition[1] + 1])
                except IndexError:
                    wR = 0
                persistence = App.persistence(persistence, orientation)
                mat[agentPosition[0]][agentPosition[1]] = 0
                mat, agentPosition, orientation = App.moveAgent(mat, agentPosition, ((wU * persistence[0]) / ((wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) + (wR * persistence[3]))), ((wD * persistence[1]) / ((wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) + (wR * persistence[3]))), ((wL * persistence[2]) / ((wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) + (wR * persistence[3]))), ((wR * persistence[3]) / ((wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) + (wR * persistence[3]))), wU, wD, wL, wR, persistence, orientation)
                # mem = App.memoryReduction(mem)
                mem[agentPosition[0]][agentPosition[1]] = mem[agentPosition[0]][agentPosition[1]] + 1
                if steps == times[0] - 1:
                    values1.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[1] - 1:
                    values2.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[2] - 1:
                    values3.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[3] - 1:
                    values4.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[4] - 1:
                    values5.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[5] - 1:
                    values6.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[6] - 1:
                    values7.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[7] - 1:
                    values8.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[8] - 1:
                    values9.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[9] - 1:
                    values10.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[10] - 1:
                    values11.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[11] - 1:
                    values12.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[12] - 1:
                    values13.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[13] - 1:
                    values14.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[14] - 1:
                    values15.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[15] - 1:
                    values16.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[16] - 1:
                    values17.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[17] - 1:
                    values18.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[18] - 1:
                    values19.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
                elif steps == times[19] - 1:
                    values20.append((agentPosition[0] - int(matrixSize / 2)) ** 2 + (agentPosition[1] - int(matrixSize / 2)) ** 2)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values1[k]
            mstd += values1[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[0], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values2[k]
            mstd += values2[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[1], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values3[k]
            mstd += values3[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[2], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values4[k]
            mstd += values4[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[3], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values5[k]
            mstd += values5[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[4], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values6[k]
            mstd += values6[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[5], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values7[k]
            mstd += values7[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[6], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values8[k]
            mstd += values8[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[7], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values9[k]
            mstd += values9[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[8], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values10[k]
            mstd += values10[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[9], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values11[k]
            mstd += values11[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[10], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values12[k]
            mstd += values12[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[11], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values13[k]
            mstd += values13[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[12], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values14[k]
            mstd += values14[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[13], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values15[k]
            mstd += values15[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[14], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values16[k]
            mstd += values16[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[15], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values17[k]
            mstd += values17[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[16], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values18[k]
            mstd += values18[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[17], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values19[k]
            mstd += values19[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[18], msd, merr)
        msd = 0.0
        mstd = 0.0
        for k in range(iterations):
            msd += values20[k]
            mstd += values20[k] ** 2
        msd = msd / iterations
        mstd = mstd / iterations
        merr = np.sqrt((mstd - msd * msd) / iterations)
        print(times[19], msd, merr)


if __name__ == "__main__":
    App.main()
