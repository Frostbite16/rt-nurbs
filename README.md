<pre align="center">

██████╗ ████████╗      ███╗   ██╗██╗   ██╗██████╗ ██████╗ ███████╗
██╔══██╗╚══██╔══╝      ████╗  ██║██║   ██║██╔══██╗██╔══██╗██╔════╝
██████╔╝   ██║   █████╗██╔██╗ ██║██║   ██║██████╔╝██████╔╝███████╗
██╔══██╗   ██║   ╚════╝██║╚██╗██║██║   ██║██╔══██╗██╔══██╗╚════██║
██║  ██║   ██║         ██║ ╚████║╚██████╔╝██║  ██║██████╔╝███████║
╚═╝  ╚═╝   ╚═╝         ╚═╝  ╚═══╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚══════╝

</pre>

> *"Offline rendering is the pursuit of optical truth, unrestricted by the illusion of time. Where parametric geometry meets mass parallel compute."*

A high-performance, offline ray tracer focused on the mathematical evaluation of NURBS (Non-Uniform Rational B-Splines) surfaces using CUDA. 

This project explores physically based rendering concepts and Data-Oriented Design (DOP) to optimize memory coalescing and mitigate thread divergence on SIMT architectures, pushing the limits of analytical surface intersection.

## Tech Stack
* **Language:** C++ / CUDA
* **Build System:** CMake
* **Windowing/Visualization:** Raylib

## Project Structure
* `src/` - Core implementation files (`.cpp`, `.cu`).
* `include/` - Header files (`.h`, `.cuh`) containing Data-Oriented structs and architecture definitions.
* `build/` - Generated build artifacts and binaries (Ignored by version control).

## Getting Started

### Prerequisites
* CMake (v3.18 or higher)
* NVIDIA CUDA Toolkit (nvcc)
* A C++ compiler compatible with your OS (GCC/Clang/MSVC)

### Building the Project
This project uses CMake to handle dependencies (like Raylib) and compilation. To build from scratch, run the following commands in the root directory:

1. Generate the build files:
   ```bash
   cmake -S . -B build


