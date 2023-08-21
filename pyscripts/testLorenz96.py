from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

# These are our constants
N = 5  # Number of variables
F = 8  # Forcing


def L96(x, t):
    """Lorenz 96 model with constant forcing"""
    # Setting up vector
    d = np.zeros(N)
    # Loops over indices (with operations and Python underflow indexing handling edge cases)
    for i in range(N):
        d[i] = (x[(i + 1) % N] - x[i - 2]) * x[i - 1] - x[i] + F
    return d


# def L96(state, t):
#     J = len(state)
#     k = np.zeros(J)

#     k[0] = (state[1] - state[J - 2]) * state[J - 1] - state[0]
#     k[1] = (state[2] - state[J - 1]) * state[0] - state[1]
#     k[J - 1] = (state[0] - state[J - 3]) * state[J - 2] - state[J - 1]

#     for j in range(2, J - 1):
#         k[j] = (state[j + 1] - state[j - 2]) * state[j - 1] - state[j]

#     return k + F


x0 = F * np.ones(N)  # Initial state (equilibrium)
x0[0] += 0.01  # Add small perturbation to the first variable
t = np.arange(0.0, 10, 0.01)

x = odeint(L96, x0, t)
print(t)
# print(x[0:3, :])
print(x)

# Plot the first three variables
fig = plt.figure()
ax = fig.add_subplot(projection="3d")
ax.plot(x[:, 0], x[:, 1], x[:, 2])
ax.set_xlabel("$x_1$")
ax.set_ylabel("$x_2$")
ax.set_zlabel("$x_3$")
plt.show()
