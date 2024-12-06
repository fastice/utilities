#
#
import strip as u
import matplotlib.pyplot as plt

def main():
    a = u.strip('/home/ian/test.dat')
    print("data")
    print(a)
    plt.plot(a[:,0],a[:,2])
    plt.show()

main()

