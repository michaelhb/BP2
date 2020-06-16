import sys, csv
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np

def main():
    reader = csv.reader(sys.stdin, delimiter = ',', quoting=csv.QUOTE_NONNUMERIC)
    data = defaultdict(lambda: defaultdict(list))
    
    for row in reader:
        tag, delta, action, time = row
        data[tag]["delta"].append(delta)
        data[tag]["action"].append(action)
        data[tag]["time"].append(time)

    ax_action = plt.subplot(211)
    ax_action.set(title=sys.argv[1], yscale='log')
    ax_action.set_xlabel("delta")
    ax_action.set_ylabel("action")
    
    ax_time = plt.subplot(212, sharex=ax_action)
    ax_time.set_ylabel("time (s)")
    
    for tag, cols in data.items():
        ax_action.plot(cols["delta"], cols["action"], label=tag)
        ax_time.plot(cols["delta"], cols["time"], label=tag)

    ax_action.legend()
    plt.show()

if __name__ == "__main__":
    main()