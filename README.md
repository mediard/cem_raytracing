
# cem
Computational Electromagnetics

## raytracing
First release is pre-alpha. Includes three files: indoor.f90, raytracing.f90 and all_data.f90. myroommm.go3 is a text file defining the floor plan with multilayered walls.

These are from about 10 years ago when I was not professionally involved with software and programming. So commenting and documentation is poor. But the algorithm and implementation are excellent. Although commenting is poor, someone familiar with ray-tracing should still be able to follow the algorithm in the raytracing module.

Here's what the tracing code does.
1) The floor plan is defined in myroommm.go3 with multiplayered walls. The complex reflection coefficients of the incident electromagnetic waves are computed based on angle of incident.
2) Image-based algorithm is used in the ray-tracing code.
3) The .TBL files are field strength.

Compile: gfortran all_data.f90 raytracing.f90 indoor.f90 -o indoor_program -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk
The -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk was necessary for me but may not be necessary for you.
Then run ./indoor_program

Soon a well-document C++ version will be released. The theory will be included. A YouTube video will be added too.
