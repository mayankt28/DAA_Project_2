from functools import cmp_to_key
import time

# Initialize the center point as (0, 0)
center = [0, 0]

# Function to calculate the orientation of three points
def calculate_orientation(a, b, c):
    result = (b[1] - a[1]) * (c[0] - b[0]) - (c[1] - b[1]) * (b[0] - a[0])
    if result == 0:
        return 0
    if result > 0:
        return 1
    return -1

# Function to merge two convex hulls
def merge_convex_hulls(a, b):
    n1, n2 = len(a), len(b)
    ia, ib = 0, 0

    # Find the rightmost point in the left convex hull
    for i in range(1, n1):
        if a[i][0] > a[ia][0]:
            ia = i

    # Find the leftmost point in the right convex hull
    for i in range(1, n2):
        if b[i][0] < b[ib][0]:
            ib = i

    inda, indb = ia, ib
    done = 0
    while not done:
        done = 1
        # Find the next point in the left convex hull
        while calculate_orientation(b[indb], a[inda], a[(inda + 1) % n1]) >= 0:
            inda = (inda + 1) % n1

        # Find the next point in the right convex hull
        while calculate_orientation(a[inda], b[indb], b[(n2 + indb - 1) % n2]) <= 0:
            indb = (indb - 1) % n2
            done = 0

    uppera, upperb = inda, indb
    inda, indb = ia, ib
    done = 0
    g = 0
    while not done:
        done = 1
        # Find the next point in the left convex hull
        while calculate_orientation(a[inda], b[indb], b[(indb + 1) % n2]) >= 0:
            indb = (indb + 1) % n2

        # Find the next point in the right convex hull
        while calculate_orientation(b[indb], a[inda], a[(n1 + inda - 1) % n1]) <= 0:
            inda = (inda - 1) % n1
            done = 0

    ret = []
    lowera, lowerb = inda, indb

    ind = uppera
    ret.append(a[uppera])
    while ind != lowera:
        ind = (ind + 1) % n1
        ret.append(a[ind])

    ind = lowerb
    ret.append(b[lowerb])
    while ind != upperb:
        ind = (ind + 1) % n2
        ret.append(b[ind])
    return ret

# Function to calculate the convex hull using brute force
def convex_hull_brute_force(a):

    # Function to compare points based on polar angles from the center
    def compare(p1, q1):
        p = [p1[0] - center[0], p1[1] - center[1]]
        q = [q1[0] - center[0], q1[1] - center[1]]
        one = None
        two = None

        # Determine the quadrant of the point
        if p[0] >= 0 and p[1] >= 0:
            one = 1
        elif p[0] <= 0 and p[1] >= 0:
            one = 2
        elif p[0] <= 0 and p[1] <= 0:
            one = 3
        else:
            one = 4

        if q[0] >= 0 and q[1] >= 0:
            two = 1
        elif q[0] <= 0 and q[1] >= 0:
            two = 2
        elif q[0] <= 0 and q[1] <= 0:
            two = 3
        else:
            two = 4

        if one != two:
            if one < two:
                return -1
            return 1
        if p[1] * q[0] < q[1] * p[0]:
            return -1
        return 1

    global center
    unique_points = set()
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            x1, x2 = a[i][0], a[j][0]
            y1, y2 = a[i][1], a[j][1]
            a1, b1, c1 = y1 - y2, x2 - x1, x1 * y2 - y1 * x2
            pos, neg = 0, 0
            for k in range(len(a)):
                if (k == i) or (k == j) or (a1 * a[k][0] + b1 * a[k][1] + c1 <= 0):
                    neg += 1
                if (k == i) or (k == j) or (a1 * a[k][0] + b1 * a[k][1] + c1 >= 0):
                    pos += 1
            if pos == len(a) or neg == len(a):
                unique_points.add(tuple(a[i]))
                unique_points.add(tuple(a[j]))

    hull_points = []
    for x in unique_points:
        hull_points.append(list(x))

    # Sorting the points in counterclockwise order
    center = [0, 0]
    n = len(hull_points)
    for i in range(n):
        center[0] += hull_points[i][0]
        center[1] += hull_points[i][1]
        hull_points[i][0] *= n
        hull_points[i][1] *= n
    hull_points = sorted(hull_points, key=cmp_to_key(compare))
    for i in range(n):
        hull_points[i] = [hull_points[i][0] / n, hull_points[i][1] / n]
    return hull_points

# Function to calculate the convex hull using the divide-and-conquer algorithm
def divide_and_conquer_convex_hull(a):

    if len(a) <= 5:
        return convex_hull_brute_force(a)

    left, right = [], []
    start = int(len(a) / 2)
    for i in range(start):
        left.append(a[i])
    for i in range(start, len(a)):
        right.append(a[i])

    left_hull = divide_and_conquer_convex_hull(left)
    right_hull = divide_and_conquer_convex_hull(right)

    return merge_convex_hulls(left_hull, right_hull)

if __name__ == '__main__':
    # Input points
    input_points = [
    (79.62, 7.59), (20.69, 29.22), (13.52, 35.68), (92.11, 95.32), (18.72, 48.98),
    (54.09, 73.68), (86.63, 80.72), (32.03, 27.01), (26.83, 41.88), (1.29, 4.42),
    (34.37, 89.09), (72.86, 1.19), (40.37, 77.53), (42.56, 24.42), (61.47, 72.09),
    (55.23, 18.01), (70.48, 54.11), (38.72, 61.64), (53.19, 52.44), (81.22, 13.79),
    (68.36, 65.91), (96.49, 43.64), (47.98, 56.55), (98.76, 88.03), (92.31, 61.59),
    (11.22, 81.11), (24.76, 38.98), (95.47, 94.89), (88.32, 32.35), (63.85, 74.22)
]

    n = len(input_points)
    input_points.sort()
    start_time = time.perf_counter_ns()
    convex_hull = divide_and_conquer_convex_hull(input_points)
    end_time = time.perf_counter_ns()

    # Calculate the elapsed time
    elapsed_time = end_time - start_time

    print(f"Elapsed Time: {elapsed_time} nano seconds")

    print('Convex Hull:')
    for point in convex_hull:
        print(int(point[0]), int(point[1]))
