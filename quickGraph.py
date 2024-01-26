from matplotlib import pyplot as plt


x= [14,
    30,
    39,
    47,
    99,
    118,
    280,
    817]

y = [0.501,
     0.411,
     0.418,
     0.506,
     0.483,
     0.481,
     0.321,
     0.195]

plt.plot(x, y, marker='o', linestyle='-')
plt.title("Average p-val vs no. communities")
plt.xlabel("No. Communities")
plt.ylabel("Average p-value")
plt.imshow()
plt.savefig("CommuntiespValue.png")
plt.show()


print("Done")