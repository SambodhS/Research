import random
import matplotlib.pyplot as plt
import numpy as np

matrixSize = 601
iterations = 5
omega1 = 0
omega2 = 0.5
epsilon = 0.025
numberOfAgents = 8
num = numberOfAgents / 2
k = 0
densitySize = 1
baselineAgents = []
testedAgents = []

times = np.logspace(0.1, 4.0, num=20)
for i in range(len(times)):
    times[i] = int(times[i])
print(times)

bvalues1 = []
bvalues2 = []
bvalues3 = []
bvalues4 = []
bvalues5 = []
bvalues6 = []
bvalues7 = []
bvalues8 = []
bvalues9 = []
bvalues10 = []
bvalues11 = []
bvalues12 = []
bvalues13 = []
bvalues14 = []
bvalues15 = []
bvalues16 = []
bvalues17 = []
bvalues18 = []
bvalues19 = []
bvalues20 = []

tvalues1 = []
tvalues2 = []
tvalues3 = []
tvalues4 = []
tvalues5 = []
tvalues6 = []
tvalues7 = []
tvalues8 = []
tvalues9 = []
tvalues10 = []
tvalues11 = []
tvalues12 = []
tvalues13 = []
tvalues14 = []
tvalues15 = []
tvalues16 = []
tvalues17 = []
tvalues18 = []
tvalues19 = []
tvalues20 = []


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

    def placeAgent(mat, mem):
        agentArray = []
        for place in range(numberOfAgents):
            done = False
            while not done:
                randomX = random.randint(int(matrixSize / 2) - densitySize, int(matrixSize / 2) + densitySize)
                randomY = random.randint(int(matrixSize / 2) - densitySize, int(matrixSize / 2) + densitySize)
                if mat[randomX][randomY] == 0:
                    if place % 2 == 0:
                        mat[randomX][randomY] = 1
                        mem[randomX][randomY] = 1
                        agentArray.append(Agent([randomX, randomY], [randomX, randomY], Main.getOrientation(), [0] * 4, 1, omega1))
                    else:
                        mat[randomX][randomY] = 2
                        mem[randomX][randomY] = 1
                        agentArray.append(Agent([randomX, randomY], [randomX, randomY], Main.getOrientation(), [0] * 4, 2, omega2))
                    done = True
        return mat, mem, agentArray


def main():
    for iteration in range(iterations):
        mat, mem, agentArray = Main.generateMatrix()
        print(iteration)
        for steps in range(int(times[19])):
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
            mem = Main.memoryReduction(mem)
            bmsd = 0
            tmsd = 0
            if steps == times[0] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues1.append(bmsd / num)
                tvalues1.append(tmsd / num)
            elif steps == times[1] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues2.append(bmsd / num)
                tvalues2.append(tmsd / num)
            elif steps == times[2] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues3.append(bmsd / num)
                tvalues3.append(tmsd / num)
            elif steps == times[3] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues4.append(bmsd / num)
                tvalues4.append(tmsd / num)
            elif steps == times[4] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues5.append(bmsd / num)
                tvalues5.append(tmsd / num)
            elif steps == times[5] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues6.append(bmsd / num)
                tvalues6.append(tmsd / num)
            elif steps == times[6] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues7.append(bmsd / num)
                tvalues7.append(tmsd / num)
            elif steps == times[7] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues8.append(bmsd / num)
                tvalues8.append(tmsd / num)
            elif steps == times[8] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues9.append(bmsd / num)
                tvalues9.append(tmsd / num)
            elif steps == times[9] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues10.append(bmsd / num)
                tvalues10.append(tmsd / num)
            elif steps == times[10] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues11.append(bmsd / num)
                tvalues11.append(tmsd / num)
            elif steps == times[11] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues12.append(bmsd / num)
                tvalues12.append(tmsd / num)
            elif steps == times[12] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues13.append(bmsd / num)
                tvalues13.append(tmsd / num)
            elif steps == times[13] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues14.append(bmsd / num)
                tvalues14.append(tmsd / num)
            elif steps == times[14] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues15.append(bmsd / num)
                tvalues15.append(tmsd / num)
            elif steps == times[15] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues16.append(bmsd / num)
                tvalues16.append(tmsd / num)
            elif steps == times[16] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues17.append(bmsd / num)
                tvalues17.append(tmsd / num)
            elif steps == times[17] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues18.append(bmsd / num)
                tvalues18.append(tmsd / num)
            elif steps == times[18] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues19.append(bmsd / num)
                tvalues19.append(tmsd / num)
            elif steps == times[19] - 1:
                for agent in agentArray:
                    if agent.designation == 1:
                        bmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                    else:
                        tmsd += (agent.agentPosition[0] - agent.startingPosition[0]) ** 2 + (agent.agentPosition[1] - agent.startingPosition[1]) ** 2
                bvalues20.append(bmsd / num)
                tvalues20.append(tmsd / num)
        # plt.imshow(mat)
        # plt.colorbar()
        # plt.show()


main()
f = open('output.txt', 'a')
msd = 0
mstd = 0
for i in range(len(bvalues1)):
    msd += bvalues1[i]
    mstd += bvalues1[i] ** 2
msd /= len(bvalues1)
mstd /= len(bvalues1)
merr = np.sqrt((mstd - msd * msd) / len(bvalues1))
f.write(str(times[0]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues2)):
    msd += bvalues2[i]
    mstd += bvalues2[i] ** 2
msd /= iterations
mstd /= iterations
merr = np.sqrt((mstd - msd * msd) / len(bvalues2))
f.write(str(times[1]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues3)):
    msd += bvalues3[i]
    mstd += bvalues3[i] ** 2
msd /= len(bvalues3)
mstd /= len(bvalues3)
merr = np.sqrt((mstd - msd * msd) / len(bvalues3))
f.write(str(times[2]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues4)):
    msd += bvalues4[i]
    mstd += bvalues4[i] ** 2
msd /= len(bvalues4)
mstd /= len(bvalues4)
merr = np.sqrt((mstd - msd * msd) / len(bvalues4))
f.write(str(times[3]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues5)):
    msd += bvalues5[i]
    mstd += bvalues5[i] ** 2
msd /= len(bvalues5)
mstd /= len(bvalues5)
merr = np.sqrt((mstd - msd * msd) / len(bvalues5))
f.write(str(times[4]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues6)):
    msd += bvalues6[i]
    mstd += bvalues6[i] ** 2
msd /= len(bvalues6)
mstd /= len(bvalues6)
merr = np.sqrt((mstd - msd * msd) / len(bvalues6))
f.write(str(times[5]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues7)):
    msd += bvalues7[i]
    mstd += bvalues7[i] ** 2
msd /= len(bvalues7)
mstd /= len(bvalues7)
merr = np.sqrt((mstd - msd * msd) / len(bvalues7))
f.write(str(times[6]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues8)):
    msd += bvalues8[i]
    mstd += bvalues8[i] ** 2
msd /= len(bvalues8)
mstd /= len(bvalues8)
merr = np.sqrt((mstd - msd * msd) / len(bvalues8))
f.write(str(times[7]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues9)):
    msd += bvalues9[i]
    mstd += bvalues9[i] ** 2
msd /= len(bvalues9)
mstd /= len(bvalues9)
merr = np.sqrt((mstd - msd * msd) / len(bvalues9))
f.write(str(times[8]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues10)):
    msd += bvalues10[i]
    mstd += bvalues10[i] ** 2
msd /= len(bvalues10)
mstd /= len(bvalues10)
merr = np.sqrt((mstd - msd * msd) / len(bvalues10))
f.write(str(times[9]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues11)):
    msd += bvalues11[i]
    mstd += bvalues11[i] ** 2
msd /= len(bvalues11)
mstd /= len(bvalues11)
merr = np.sqrt((mstd - msd * msd) / len(bvalues11))
f.write(str(times[10]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues12)):
    msd += bvalues12[i]
    mstd += bvalues12[i] ** 2
msd /= len(bvalues12)
mstd /= len(bvalues12)
merr = np.sqrt((mstd - msd * msd) / len(bvalues12))
f.write(str(times[11]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues13)):
    msd += bvalues13[i]
    mstd += bvalues13[i] ** 2
msd /= len(bvalues13)
mstd /= len(bvalues13)
merr = np.sqrt((mstd - msd * msd) / len(bvalues13))
f.write(str(times[12]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues14)):
    msd += bvalues14[i]
    mstd += bvalues14[i] ** 2
msd /= len(bvalues14)
mstd /= len(bvalues14)
merr = np.sqrt((mstd - msd * msd) / len(bvalues14))
f.write(str(times[13]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues15)):
    msd += bvalues15[i]
    mstd += bvalues15[i] ** 2
msd /= len(bvalues15)
mstd /= len(bvalues15)
merr = np.sqrt((mstd - msd * msd) / len(bvalues15))
f.write(str(times[14]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues16)):
    msd += bvalues16[i]
    mstd += bvalues16[i] ** 2
msd /= len(bvalues16)
mstd /= len(bvalues16)
merr = np.sqrt((mstd - msd * msd) / len(bvalues16))
f.write(str(times[15]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues17)):
    msd += bvalues17[i]
    mstd += bvalues17[i] ** 2
msd /= len(bvalues17)
mstd /= len(bvalues17)
merr = np.sqrt((mstd - msd * msd) / len(bvalues17))
f.write(str(times[16]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues18)):
    msd += bvalues18[i]
    mstd += bvalues18[i] ** 2
msd /= len(bvalues18)
mstd /= len(bvalues18)
merr = np.sqrt((mstd - msd * msd) / len(bvalues18))
f.write(str(times[17]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues19)):
    msd += bvalues19[i]
    mstd += bvalues19[i] ** 2
msd /= len(bvalues19)
mstd /= len(bvalues19)
merr = np.sqrt((mstd - msd * msd) / len(bvalues19))
f.write(str(times[18]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(bvalues20)):
    msd += bvalues20[i]
    mstd += bvalues20[i] ** 2
msd /= len(bvalues20)
mstd /= len(bvalues20)
merr = np.sqrt((mstd - msd * msd) / len(bvalues20))
f.write(str(times[19]) + ", " + str(msd) + ", " + str(merr) + "\n\n")

msd = 0
mstd = 0
for i in range(len(tvalues1)):
    msd += tvalues1[i]
    mstd += tvalues1[i] ** 2
msd /= len(tvalues1)
mstd /= len(tvalues1)
merr = np.sqrt((mstd - msd * msd) / len(tvalues1))
f.write(str(times[0]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues2)):
    msd += tvalues2[i]
    mstd += tvalues2[i] ** 2
msd /= len(tvalues2)
mstd /= len(tvalues2)
merr = np.sqrt((mstd - msd * msd) / len(tvalues2))
f.write(str(times[1]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues3)):
    msd += tvalues3[i]
    mstd += tvalues3[i] ** 2
msd /= len(tvalues3)
mstd /= len(tvalues3)
merr = np.sqrt((mstd - msd * msd) / len(tvalues3))
f.write(str(times[2]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues4)):
    msd += tvalues4[i]
    mstd += tvalues4[i] ** 2
msd /= len(tvalues4)
mstd /= len(tvalues4)
merr = np.sqrt((mstd - msd * msd) / len(tvalues4))
f.write(str(times[3]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues5)):
    msd += tvalues5[i]
    mstd += tvalues5[i] ** 2
msd /= len(tvalues5)
mstd /= len(tvalues5)
merr = np.sqrt((mstd - msd * msd) / len(tvalues5))
f.write(str(times[4]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues6)):
    msd += tvalues6[i]
    mstd += tvalues6[i] ** 2
msd /= len(tvalues6)
mstd /= len(tvalues6)
merr = np.sqrt((mstd - msd * msd) / len(tvalues6))
f.write(str(times[5]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues7)):
    msd += tvalues7[i]
    mstd += tvalues7[i] ** 2
msd /= len(tvalues7)
mstd /= len(tvalues7)
merr = np.sqrt((mstd - msd * msd) / len(tvalues7))
f.write(str(times[6]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues8)):
    msd += tvalues8[i]
    mstd += tvalues8[i] ** 2
msd /= len(tvalues8)
mstd /= len(tvalues8)
merr = np.sqrt((mstd - msd * msd) / len(tvalues8))
f.write(str(times[7]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues9)):
    msd += tvalues9[i]
    mstd += tvalues9[i] ** 2
msd /= len(tvalues9)
mstd /= len(tvalues9)
merr = np.sqrt((mstd - msd * msd) / len(tvalues9))
f.write(str(times[8]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues10)):
    msd += tvalues10[i]
    mstd += tvalues10[i] ** 2
msd /= len(tvalues10)
mstd /= len(tvalues10)
merr = np.sqrt((mstd - msd * msd) / len(tvalues10))
f.write(str(times[9]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues11)):
    msd += tvalues11[i]
    mstd += tvalues11[i] ** 2
msd /= len(tvalues11)
mstd /= len(tvalues11)
merr = np.sqrt((mstd - msd * msd) / len(tvalues11))
f.write(str(times[10]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues12)):
    msd += tvalues12[i]
    mstd += tvalues12[i] ** 2
msd /= len(tvalues12)
mstd /= len(tvalues12)
merr = np.sqrt((mstd - msd * msd) / len(tvalues12))
f.write(str(times[11]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues13)):
    msd += tvalues13[i]
    mstd += tvalues13[i] ** 2
msd /= len(tvalues13)
mstd /= len(tvalues13)
merr = np.sqrt((mstd - msd * msd) / len(tvalues13))
f.write(str(times[12]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues14)):
    msd += tvalues14[i]
    mstd += tvalues14[i] ** 2
msd /= len(tvalues14)
mstd /= len(tvalues14)
merr = np.sqrt((mstd - msd * msd) / len(tvalues14))
f.write(str(times[13]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues15)):
    msd += tvalues15[i]
    mstd += tvalues15[i] ** 2
msd /= len(tvalues15)
mstd /= len(tvalues15)
merr = np.sqrt((mstd - msd * msd) / len(tvalues15))
f.write(str(times[14]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues16)):
    msd += tvalues16[i]
    mstd += tvalues16[i] ** 2
msd /= len(tvalues16)
mstd /= len(tvalues16)
merr = np.sqrt((mstd - msd * msd) / len(tvalues16))
f.write(str(times[15]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues17)):
    msd += tvalues17[i]
    mstd += tvalues17[i] ** 2
msd /= len(tvalues17)
mstd /= len(tvalues17)
merr = np.sqrt((mstd - msd * msd) / len(tvalues17))
f.write(str(times[16]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues18)):
    msd += tvalues18[i]
    mstd += tvalues18[i] ** 2
msd /= len(tvalues18)
mstd /= len(tvalues18)
merr = np.sqrt((mstd - msd * msd) / len(tvalues18))
f.write(str(times[17]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues19)):
    msd += tvalues19[i]
    mstd += tvalues19[i] ** 2
msd /= len(tvalues19)
mstd /= len(tvalues19)
merr = np.sqrt((mstd - msd * msd) / len(tvalues19))
f.write(str(times[18]) + ", " + str(msd) + ", " + str(merr) + "\n")
msd = 0
mstd = 0
for i in range(len(tvalues20)):
    msd += tvalues20[i]
    mstd += tvalues20[i] ** 2
msd /= len(tvalues20)
mstd /= len(tvalues20)
merr = np.sqrt((mstd - msd * msd) / len(tvalues20))
f.write(str(times[19]) + ", " + str(msd) + ", " + str(merr) + "\n")
