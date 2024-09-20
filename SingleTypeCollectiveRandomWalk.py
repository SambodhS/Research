import random
import numpy as np

matrixSize = 501
iterations = 5
epsilon = 0.025
numberOfAgents = 8
densitySize = 1
k = 0

# 0.3, 0.05, 0.5omega for emergent collective motion

times = np.logspace(0.1, 4.0, num=20)
for i in range(len(times)):
    times[i] = int(times[i])

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


class Agent:
    def __init__(self, agentPosition, startingPosition, orientation, persistence, designation, omega):
        self.agentPosition = agentPosition
        self.startingPosition = startingPosition
        self.orientation = orientation
        self.persistence = persistence
        self.designation = designation
        self.omega = omega

    def moveAgent(self, mat, mem, pU, pD, pL, pR, wU, wD, wL, wR):
        r = random.uniform(0, 1)
        mat[self.agentPosition[0]][self.agentPosition[1]] = 0
        if self.agentPosition[0] == 0 or self.agentPosition[0] == matrixSize - 1 or self.agentPosition[1] == 0 or self.agentPosition == matrixSize - 1:
            print("Oh no!")
        if r <= pU:
            if self.agentPosition[0] != 0 and mat[self.agentPosition[0] - 1][self.agentPosition[1]] == 0:
                mat[self.agentPosition[0] - 1][self.agentPosition[1]] = self.designation
                mem[self.agentPosition[0] - 1][self.agentPosition[1]] += 1
                self.agentPosition[0] = self.agentPosition[0] - 1
                self.orientation = "U"
            elif self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0:
                if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, 0, (wD * self.persistence[1]) / (wD * self.persistence[1] + wL * self.persistence[2] + wR * self.persistence[3]),
                                       (wL * self.persistence[2]) / (wD * self.persistence[1] + wL * self.persistence[2] + wR * self.persistence[3]),
                                       (wR * self.persistence[3]) / (wD * self.persistence[1] + wL * self.persistence[2] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 0, (wD * self.persistence[1]) / (wD * self.persistence[1] + wL * self.persistence[2]),
                                       (wL * self.persistence[2]) / (wD * self.persistence[1] + wL * self.persistence[2]), 0, wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, 0, (wD * self.persistence[1]) / (wD * self.persistence[1] + wR * self.persistence[3]), 0,
                                       (wR * self.persistence[3]) / (wD * self.persistence[1] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 0, 1, 0, 0, wU, wD, wL, wR)
            else:
                if (self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0) and (self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0):
                    self.moveAgent(mat, mem, 0, 0, (wL * self.persistence[2]) / (wL * self.persistence[2] + wR * self.persistence[3]),
                                   (wR * self.persistence[3]) / (wL * self.persistence[2] + wR * self.persistence[3]), wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                        self.moveAgent(mat, mem, 0, 0, 1, 0, wU, wD, wL, wR)
                    else:
                        if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                            self.moveAgent(mat, mem, 0, 0, 0, 1, wU, wD, wL, wR)
        elif pU < r <= pU + pD:
            if self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0:
                mat[self.agentPosition[0] + 1][self.agentPosition[1]] = self.designation
                mem[self.agentPosition[0] + 1][self.agentPosition[1]] += 1
                self.agentPosition[0] = self.agentPosition[0] + 1
                self.orientation = "D"
            elif self.agentPosition[0] != 0 and mat[self.agentPosition[0] - 1][self.agentPosition[1]] == 0:
                if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wL * self.persistence[2] + wR * self.persistence[3]),
                                       0, (wL * self.persistence[2]) / (wU * self.persistence[0] + wL * self.persistence[2] + wR * self.persistence[3]),
                                       (wR * self.persistence[3]) / (wU * self.persistence[0] + wL * self.persistence[2] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wL * self.persistence[2]), 0,
                                       (wL * self.persistence[2]) / (wU * self.persistence[0] + wL * self.persistence[2]), 0, wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wR * self.persistence[3]), 0, 0,
                                       (wR * self.persistence[3]) / (wU * self.persistence[0] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 1, 0, 0, 0, wU, wD, wL, wR)
            else:
                if (self.agentPosition != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0) and (self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0):
                    self.moveAgent(mat, mem, 0, 0, (wL * self.persistence[2]) / (wL * self.persistence[2] + wR * self.persistence[3]),
                                   (wR * self.persistence[3]) / (wL * self.persistence[2] + wR * self.persistence[3]), wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                        self.moveAgent(mat, mem, 0, 0, 1, 0, wU, wD, wL, wR)
                    else:
                        if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                            self.moveAgent(mat, mem, 0, 0, 0, 1, wU, wD, wL, wR)
        elif pU + pD < r <= pU + pD + pL:
            if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                mat[self.agentPosition[0]][self.agentPosition[1] - 1] = self.designation
                mem[self.agentPosition[0]][self.agentPosition[1] - 1] += 1
                self.agentPosition[1] = self.agentPosition[1] - 1
                self.orientation = "L"
            elif self.agentPosition[0] != 0 and mat[self.agentPosition[0] - 1][self.agentPosition[1]] == 0:
                if self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wD * self.persistence[1] + wR * self.persistence[3]),
                                       (wD * self.persistence[1]) / (wU * self.persistence[0] + wD * self.persistence[1] + wR * self.persistence[3]), 0,
                                       (wR * self.persistence[3]) / (wU * self.persistence[0] + wD * self.persistence[1] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wD * self.persistence[1]),
                                  (wD * self.persistence[1]) / (wU * self.persistence[0] + wD * self.persistence[1]), 0, 0, wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wR * self.persistence[3]), 0, 0,
                                       (wR * self.persistence[3]) / (wU * self.persistence[0] + wR * self.persistence[3]), wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 1, 0, 0, 0, wU, wD, wL, wR)
            else:
                if (self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0) and (self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0):
                    self.moveAgent(mat, mem, 0, (wD * self.persistence[1]) / (wD * self.persistence[1] + wR * self.persistence[3]), 0,
                                   (wR * self.persistence[3]) / (wD * self.persistence[1] + wR * self.persistence[3]), wU, wD, wL, wR)
                else:
                    if self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0:
                        self.moveAgent(mat, mem, 0, 1, 0, 0, wU, wD, wL, wR)
                    else:
                        if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                            self.moveAgent(mat, mem, 0, 0, 0, 1, wU, wD, wL, wR)
        elif pU + pD + pL < r <= pU + pD + pL + pR:
            if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                mat[self.agentPosition[0]][self.agentPosition[1] + 1] = self.designation
                mem[self.agentPosition[0]][self.agentPosition[1] + 1] += 1
                self.agentPosition[1] = self.agentPosition[1] + 1
                self.orientation = "R"
            elif self.agentPosition[0] != 0 and mat[self.agentPosition[0] - 1][self.agentPosition[1]] == 0:
                if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                    if self.agentPosition[1] != matrixSize - 1 and mat[self.agentPosition[0]][self.agentPosition[1] + 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wD * self.persistence[1] + wL * self.persistence[2]),
                                       (wD * self.persistence[1]) / (wU * self.persistence[0] + wD * self.persistence[1] + wL * self.persistence[2]),
                                       (wL * self.persistence[2]) / (wU * self.persistence[0] + wD * self.persistence[1] + wL * self.persistence[2]), 0, wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wD * self.persistence[1]),
                                       (wD * self.persistence[1]) / (wU * self.persistence[0] + wD * self.persistence[1]), 0, 0, wU, wD, wL, wR)
                else:
                    if self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0:
                        self.moveAgent(mat, mem, (wU * self.persistence[0]) / (wU * self.persistence[0] + wL * self.persistence[2]), 0,
                                       (wL * self.persistence[2]) / (wU * self.persistence[0] + wL * self.persistence[2]), 0, wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 1, 0, 0, 0, wU, wD, wL, wR)
            else:
                if (self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0) and (self.agentPosition[1] != 0 and mat[self.agentPosition[0]][self.agentPosition[1] - 1] == 0):
                    self.moveAgent(mat, mem, 0, (wD * self.persistence[1]) / (wD * self.persistence[1] + wL * self.persistence[2]),
                                   (wL * self.persistence[2]) / (wD * self.persistence[1] + wL * self.persistence[2]), 0, wU, wD, wL, wR)
                else:
                    if self.agentPosition[0] != matrixSize - 1 and mat[self.agentPosition[0] + 1][self.agentPosition[1]] == 0:
                        self.moveAgent(mat, mem, 0, 1, 0, 0, wU, wD, wL, wR)
                    else:
                        self.moveAgent(mat, mem, 0, 0, 1, 0, wU, wD, wL, wR)
        return mat, mem, self.agentPosition

    def calculatePersistence(self):
        if self.orientation == "U":
            self.persistence[0] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif self.orientation == "D":
            self.persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[1] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif self.orientation == "L":
            self.persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[2] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[3] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        elif self.orientation == "R":
            self.persistence[0] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[1] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[2] = (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
            self.persistence[3] = 1 - 3 * (np.exp(-k) / (3 * np.exp(-k) + np.exp(k)))
        return self.persistence

    @staticmethod
    def weights(strength, omega):
        V = omega * strength
        if V > 709:
            V = 709
        weight = np.exp(V)
        return weight


class Main:
    @staticmethod
    def generateMatrix():
        mat = [[0] * matrixSize for _ in range(matrixSize)]
        mem = [[0] * matrixSize for _ in range(matrixSize)]
        mat, mem, agentArray = Main.placeAgent(mat, mem)
        return mat, mem, agentArray

    def print2D(mat):
        for row in mat:
            print(row)

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

    def memoryReduction(mem):
        i = 0
        while i < len(mem):
            j = 0
            while j < len(mem[i]):
                if mem[i][j] > 0:
                    mem[i][j] -= epsilon * mem[i][j]
                j += 1
            i += 1
        return mem

    def placeAgent(mat, mem):
        agentArray = []
        for place in range(numberOfAgents):
            done = False
            while not done:
                randomX = random.randint(int(matrixSize / 2) - densitySize, int(matrixSize / 2) + densitySize)
                randomY = random.randint(int(matrixSize / 2) - densitySize, int(matrixSize / 2) + densitySize)
                if mat[randomX][randomY] == 0:
                    mat[randomX][randomY] = 1
                    agentArray.append(Agent([randomX, randomY], [randomX, randomY], Main.getOrientation(), [0] * 4, 1, 0.5))
                    done = True
        return mat, mem, agentArray


def main():
    for iteration in range(iterations):
        print(iteration)
        mat, mem, agentArray = Main.generateMatrix()
        for step in range(int(times[19])):
            for agent in agentArray:
                try:
                    wU = agent.weights(mem[agent.agentPosition[0] - 1][agent.agentPosition[1]], agent.omega)
                except IndexError:
                    wU = 0
                try:
                    wD = agent.weights(mem[agent.agentPosition[0] + 1][agent.agentPosition[1]], agent.omega)
                except IndexError:
                    wD = 0
                try:
                    wL = agent.weights(mem[agent.agentPosition[0]][agent.agentPosition[1] - 1], agent.omega)
                except IndexError:
                    wL = 0
                try:
                    wR = agent.weights(mem[agent.agentPosition[0]][agent.agentPosition[1] + 1], agent.omega)
                except IndexError:
                    wR = 0
                persistence = agent.calculatePersistence()
                mat, mem, agentPosition = agent.moveAgent(mat, mem, (wU * persistence[0]) / (
                        (wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) +
                        (wR * persistence[3])), (wD * persistence[1]) / (
                        (wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) +
                        (wR * persistence[3])), (wL * persistence[2]) / (
                        (wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) +
                        (wR * persistence[3])), (wR * persistence[3]) / (
                        (wU * persistence[0]) + (wD * persistence[1]) + (wL * persistence[2]) +
                        (wR * persistence[3])), wU, wD, wL, wR)
                agent.agentPosition = agentPosition
            if step == times[0] - 1:
                for agent in agentArray:
                    values1.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[1] - 1:
                for agent in agentArray:
                    values2.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[2] - 1:
                for agent in agentArray:
                    values3.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[3] - 1:
                for agent in agentArray:
                    values4.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[4] - 1:
                for agent in agentArray:
                    values5.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[5] - 1:
                for agent in agentArray:
                    values6.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[6] - 1:
                for agent in agentArray:
                    values7.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[7] - 1:
                for agent in agentArray:
                    values8.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[8] - 1:
                for agent in agentArray:
                    values9.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[9] - 1:
                for agent in agentArray:
                    values10.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[10] - 1:
                for agent in agentArray:
                    values11.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[11] - 1:
                for agent in agentArray:
                    values12.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[12] - 1:
                for agent in agentArray:
                    values13.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[13] - 1:
                for agent in agentArray:
                    values14.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[14] - 1:
                for agent in agentArray:
                    values15.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[15] - 1:
                for agent in agentArray:
                    values16.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[16] - 1:
                for agent in agentArray:
                    values17.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[17] - 1:
                for agent in agentArray:
                    values18.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[18] - 1:
                for agent in agentArray:
                    values19.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            elif step == times[19] - 1:
                for agent in agentArray:
                    values20.append((agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2)
            mem = Main.memoryReduction(mem)
        # plt.imshow(mat)
        # plt.colorbar()
        # plt.show()


main()
f = open('output.txt', 'a')
msd = 0
mstd = 0
for i in range(len(values1)):
    msd += values1[i]
    mstd += values1[i] ** 2
msd /= len(values1)
mstd /= len(values1)
merr = np.sqrt((mstd - msd * msd) / len(values1))
f.write(str(times[0]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values2)):
    msd += values2[i]
    mstd += values2[i] ** 2
msd /= len(values2)
mstd /= len(values2)
merr = np.sqrt((mstd - msd * msd) / len(values2))
f.write(str(times[1]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values3)):
    msd += values3[i]
    mstd += values3[i] ** 2
msd /= len(values3)
mstd /= len(values3)
merr = np.sqrt((mstd - msd * msd) / len(values3))
f.write(str(times[2]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values4)):
    msd += values4[i]
    mstd += values4[i] ** 2
msd /= len(values4)
mstd /= len(values4)
merr = np.sqrt((mstd - msd * msd) / len(values4))
f.write(str(times[3]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values5)):
    msd += values5[i]
    mstd += values5[i] ** 2
msd /= len(values5)
mstd /= len(values5)
merr = np.sqrt((mstd - msd * msd) / len(values5))
f.write(str(times[4]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values6)):
    msd += values6[i]
    mstd += values6[i] ** 2
msd /= len(values6)
mstd /= len(values6)
merr = np.sqrt((mstd - msd * msd) / len(values6))
f.write(str(times[5]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values7)):
    msd += values7[i]
    mstd += values7[i] ** 2
msd /= len(values7)
mstd /= len(values7)
merr = np.sqrt((mstd - msd * msd) / len(values7))
f.write(str(times[6]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values8)):
    msd += values8[i]
    mstd += values8[i] ** 2
msd /= len(values8)
mstd /= len(values8)
merr = np.sqrt((mstd - msd * msd) / len(values8))
f.write(str(times[7]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values9)):
    msd += values9[i]
    mstd += values9[i] ** 2
msd /= len(values9)
mstd /= len(values9)
merr = np.sqrt((mstd - msd * msd) / len(values9))
f.write(str(times[8]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values10)):
    msd += values10[i]
    mstd += values10[i] ** 2
msd /= len(values10)
mstd /= len(values10)
merr = np.sqrt((mstd - msd * msd) / len(values10))
f.write(str(times[9]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values11)):
    msd += values11[i]
    mstd += values11[i] ** 2
msd /= len(values11)
mstd /= len(values11)
merr = np.sqrt((mstd - msd * msd) / len(values11))
f.write(str(times[10]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values12)):
    msd += values12[i]
    mstd += values12[i] ** 2
msd /= len(values12)
mstd /= len(values12)
merr = np.sqrt((mstd - msd * msd) / len(values12))
f.write(str(times[11]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values13)):
    msd += values13[i]
    mstd += values13[i] ** 2
msd /= len(values13)
mstd /= len(values13)
merr = np.sqrt((mstd - msd * msd) / len(values13))
f.write(str(times[12]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values14)):
    msd += values14[i]
    mstd += values14[i] ** 2
msd /= len(values14)
mstd /= len(values14)
merr = np.sqrt((mstd - msd * msd) / len(values14))
f.write(str(times[13]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values15)):
    msd += values15[i]
    mstd += values15[i] ** 2
msd /= len(values15)
mstd /= len(values15)
merr = np.sqrt((mstd - msd * msd) / len(values15))
f.write(str(times[14]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values16)):
    msd += values16[i]
    mstd += values16[i] ** 2
msd /= len(values16)
mstd /= len(values16)
merr = np.sqrt((mstd - msd * msd) / len(values16))
f.write(str(times[15]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values17)):
    msd += values17[i]
    mstd += values17[i] ** 2
msd /= len(values17)
mstd /= len(values17)
merr = np.sqrt((mstd - msd * msd) / len(values17))
f.write(str(times[16]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values18)):
    msd += values18[i]
    mstd += values18[i] ** 2
msd /= len(values18)
mstd /= len(values18)
merr = np.sqrt((mstd - msd * msd) / len(values18))
f.write(str(times[17]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values19)):
    msd += values19[i]
    mstd += values19[i] ** 2
msd /= len(values19)
mstd /= len(values19)
merr = np.sqrt((mstd - msd * msd) / len(values19))
f.write(str(times[18]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(values20)):
    msd += values20[i]
    mstd += values20[i] ** 2
msd /= len(values20)
mstd /= len(values20)
merr = np.sqrt((mstd - msd * msd) / len(values20))
f.write(str(times[19]) + ", " + str(msd) + ", " + str(merr) + "\n")
