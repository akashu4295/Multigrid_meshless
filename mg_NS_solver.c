// Author :  Akash Unnikrishnan and Prof. Surya Pratap Vanka
// Affiliation : Indian Institute of Technology Gandhinagar and University of Illinois at Urbana Champaign
// compile command: nvcc multigrid_heat_conduction.c -lm
// Run command: ./a.out

///////////////////////////////////////////////////////////////////////////////
////////////// Header files

// #include "header_files/fractional_step_solver.h"
#include "header_files/timple_solver.h"
#include "header_files/openACC_functions.h"
#include "init/init_TC.c"

///////////////////////////////////////////////////////////////////////////////
////////////// Main Program

int main()
{
    clock_t clock_start = clock(), clock_program_begin = clock();
    PointStructure* myPointStruct;
    FieldVariables *field;
    double steady_state_error;
    FILE *file1;
    FILE *file2;

///// Read Input File Name and Parameters
    read_flow_parameters("flow_parameters.csv");     // Read the parameters from the file
    read_grid_filenames(&myPointStruct, "grid_filenames.csv", &parameters.num_levels);     // Read the parameters from the file

///// Read Mesh data on all levels and compute derivative matrices
    read_complete_mesh_data(myPointStruct, parameters.num_levels);
    printf("Time taken to read the grids and flow parameters: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    for (int ii = 0; ii<parameters.num_levels ; ii = ii +1)
        create_derivative_matrices(&myPointStruct[ii]);
    printf("Time taken to create derivative matrices: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

    clock_start = clock();    // Start the clock
    if(parameters.test>0){
        test_derivatives (myPointStruct, parameters.num_levels, parameters.dimension);
        printf("Time taken to test the derivatives: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);
    }

////////////// Setting up the boudary condition (Lid driven cavity) /////////////
    
    AllocateMemoryFieldVariables(&field, myPointStruct, parameters.num_levels);

    initial_conditions(myPointStruct, field, 1);
    boundary_conditions(myPointStruct, field, 1);
    parameters.dt = calculate_dt(&myPointStruct[0]);
    //parameters.dt =  0.01;
    printf("Time to setup the problem in cpu: %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);

////////////// Copy data to GPU memory /////////////
    clock_start = clock();    // Start the clock
    copyin_essentials_to_gpu(myPointStruct);
    copyin_field_to_gpu(field, myPointStruct);
    copyin_parameters_to_gpu();

    printf("Time taken to copy data to GPU: %lf\n", (double)(clock()-clock_start)/CLOCKS_PER_SEC);

////////////// Time stepping ////////////// 
//////////// Copy data back and check the rms error /////
    
    file2 = fopen("Convergence.csv", "w"); // Write data to a file

    for (int it = 0; it<parameters.num_time_steps; it++ ) 
    {
    // if (steady_state_error < 1.0e-6) parameters.iter_timple = 5;
        // steady_state_error = fractional_step_explicit(myPointStruct, field);
       steady_state_error = time_implicit_solver(myPointStruct, field);
        printf("Time step: %d, Steady state error: %e\n", it, steady_state_error);
        
        # pragma acc update host(field[0])
        
        fprintf(file2,"%d, %e\n", it, steady_state_error);
        fflush(file2);
        
        if (steady_state_error < parameters.steady_state_tolerance){
            printf("Converged at time step: %d\n", it);
            break;
        }
        if ((it % parameters.write_interval == 0) || (it == parameters.num_time_steps-1)){
            file1 = fopen("Solution.csv", "w"); // Write data to a file
            for (int i = 0; i < myPointStruct[0].num_nodes; i++)
                fprintf(file1, "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n", myPointStruct[0].x[i], myPointStruct[0].y[i], myPointStruct[0].z[i], field[0].u[i], field[0].v[i], field[0].w[i], field[0].p_old[i]);
            fflush(file1);	
            fclose(file1);
        }
    } 
    
    printf("Time_step, dt : %lf\n",parameters.dt);
    printf("Polynomial degree: %d\n",myPointStruct[0].poly_degree);
    printf("Time for execution (total): %lf\n", (double)(clock()-clock_program_begin)/CLOCKS_PER_SEC);
    
    // Clear read data
    free_PointStructure(myPointStruct, parameters.num_levels);
    free_field(field, parameters.num_levels);

    fclose(file2); 
    return 0;
} 
