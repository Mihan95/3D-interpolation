HEADERS       += funcs.h scene3D.h solver.h solver_thread.h
                
SOURCES       = main.cpp \
                solver.cpp \
                fill_matrix.cpp \
                fill_vector.cpp \
                approximation.cpp \
                scene3D.cpp \

QT += gui core widgets opengl