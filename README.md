# tlb_exact

A Benders decomposition algorithm to solve the target-level bicycle problem for large-scale bike sharing systems. The problem was introduced in this [Technical Report](https://www.cirrelt.ca/documentstravail/cirrelt-2025-02.pdf).

## Building the code in Linux

1. Clone the repository and add executable permission to a script that will call CMake for you:

```shell
git clone https://github.com/lucasparada20/tlb_exact.git
cd sbrp
chmod u+x cmake_script_exact.sh
```
2. The compiler needs an absolute path to your installed Cplex library. To provide the path, go into src/CMakeLists.txt and edit the following line:

```cmake
set(CPLEX_DIR "/some/path/to/Cplex")
```

3. The compiler also need an absolute path to your installed [LEMON](https://lemon.cs.elte.hu/trac/lemon) library. To provide the path, go into src/CMakeLists.txt and edit the following line:

```cmake
set(LEMON_DIR "/some/path/to/lemon-1.3.1")
```
LEMON is a very efficiently library to solve network optimization problems. In the case of the TLB, it is used to solve the scenario min-cost flow subproblem.

4. Build the code by typing:

```bash
./cmake_script_exact.sh
```

The TLB uses `OpenMP` which is enabled by default in the compiler.

## Instances

The instances can be found in this [Google drive link](https://drive.google.com/file/d/1Q-0E389K-WTVqK05zVU2rITzFuYCQLoM/view?usp=sharing). These instances are very large-sized. Smaller instances can be found in my other repo for the [TLB generator](https://github.com/lucasparada20/tlb_generator).
