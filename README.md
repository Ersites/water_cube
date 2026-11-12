# Water cube

Functions to be used in a voxel water simulation.

## Convention

We assume a unit cube, between (0, 0, 0) and (1, 1, 1).

G is a normalized vector denoting the normal of the plane P. It can also be thought of as a vector opposite to the gravity vector. By being "below the plane" we mean points on the side of the P opposite to the normal, so along the gravity vector.

Max height H_max is the height of the cube when projected onto G. Its value is between 1 and sqrt(3).
The cube contains a volume V between 0 and 1 of a water (or fluid in general). For given G and V we can think of the water as settling in the cube. Its surface will eventually form the plane P, which is what we are interested in.

H is the height of the water. It is the distance of the lowest vertex of the cube to the water plane. It can be also thought of as projecting both P and the cube onto G and calculating the distance between the projection of P and the projection of the lowest point of the cube. For a given G, H will be between 0 and H_max.

We can create functions V(H) and H(V) where G is a constant.

## License
This project is licensed under the MIT License ([LICENSE](LICENSE)).