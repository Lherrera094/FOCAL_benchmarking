#include "save_data.h"

void create_folder( gridConfiguration *gridCfg, saveData *saveDCfg ){
    //{{{

    printf("-------------------Storage Folder------------------- \n");
    simulation_folder( projectPath );
    data_folder( projectPath, foldername );
    copyJSON( projectPath, foldername );

}//}}}


void simulation_folder(const char *path){
    //{{{
    
    struct stat st = {0};

    /*Checks if directory exists*/
    if( stat(path, &st) == -1){
        //Directory does not exists. Create it.
        if(mkdir(path, 0700) == 0){
            printf("Main project folder created successfully.\n");
        }else{
            printf("Error creating directory: %s\n", path);
            return;
        }
    }
}//}}}


void data_folder(const char *path, const char *folder_name){
    //{{{
    
    char fullPath[PATH_MAX];

    // Create the full directory path and check for buffer overflow
    if (snprintf(fullPath, sizeof(fullPath), "%s/%s", path, folder_name) >= sizeof(fullPath)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    struct stat st = {0};

    /*Checks if directory exists.*/
    if( stat(fullPath, &st) == -1){
        //Directory does not exists. Create it.
        if( mkdir(fullPath, 0700) == 0){
            printf("%s folder created successfully. \n", folder_name);
        }else{
            printf("Error creating directory: %s\n", folder_name);
            return;
        }
    }else{
        printf("%s already exists.\n", folder_name);
    }
}//}}}


void copyJSON(const char *path, const char *folder_name){
    //{{{
    
    char destination[PATH_MAX];

    //Read the source file
    FILE *srcFile = fopen("input_FOCAL.json", "rb");
    if (srcFile == NULL) {
        perror("Error opening source file");
        return;
    }

    // Create the full directory path and check for buffer overflow
    if (snprintf(destination, sizeof(destination), "%s/%s/input_FOCAL.json", path, folder_name) >= sizeof(destination)) {
        fprintf(stderr, "Error: Directory path is too long.\n");
        return;
    }

    //open destination file
    FILE *destFile = fopen(destination, "wb");
    if (destFile == NULL) {
        perror("Error opening destination file");
        fclose(srcFile);
        return;
    }

    char buffer[1024];
    size_t bytesRead;
    while ((bytesRead = fread(buffer, 1, sizeof(buffer), srcFile)) > 0) {
        fwrite(buffer, 1, bytesRead, destFile);
    }

    fclose(srcFile);
    fclose(destFile);

    printf("JSON file saved.\n");
}//}}}

//Functions to save data
void save_SimData(  gridConfiguration *gridCfg,
                    beamAntennaConfiguration *beamCfg, 
                    saveData *saveDCfg,
                    double n_e[NX/2][NY/2][NZ/2],
                    double J_B0[NX][NY][NZ] ){

    //Char values as directions to the correct folder
    char fullDir[PATH_MAX], filename_hdf5[PATH_MAX];

    sprintf(fullDir,"%s/%s/", projectPath, foldername);
    
    //Append the name of the files
    sprintf( filename_hdf5, "%s%s", fullDir, file_hdf5);

    /*Save data into path directory*/
    save_data_toHDF5( gridCfg, beamCfg, filename_hdf5 , n_e, J_B0 );
    
}

int save_data_toHDF5(   gridConfiguration *gridCfg,
                        beamAntennaConfiguration *beamCfg,
                        char filename_hdf5[],
                        double n_e[NX/2][NY/2][NZ/2],
                        double J_B0[NX][NY][NZ] ){
    
    int ii, jj , kk;
    char dSet_name[PATH_MAX];
    // used when writing data into hdf5-files
    double (*data2save)[NY/2][NZ/2]     = calloc(NX/2, sizeof *data2save);

    writeConfig2HDF( gridCfg, beamCfg, filename_hdf5 );

    // density
    sprintf( dSet_name, "n_e" );
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, dSet_name, n_e) ) ;

    // background magnetic field
    // B0x: even-odd-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii  ][jj+1][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0x", data2save) ) ;

    set2zero_3D( NX/2, NY/2, NZ/2, data2save );
    // B0y: odd-even-odd
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj  ][kk+1];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0y", data2save) ) ;

    set2zero_3D( NX/2, NY/2, NZ/2, data2save );
    // B0z: odd-odd-even
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
    for (ii=0 ; ii<NX ; ii+=2) {
        for (jj=0 ; jj<NY ; jj+=2) {
            for (kk=0 ; kk<NZ ; kk+=2) {
                data2save[(ii/2)][(jj/2)][(kk/2)] = J_B0[ii+1][jj+1][kk  ];
            }
        }
    }
    printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, "B0z", data2save) ) ;

    free( data2save );
    printf( "Free data2save allocated memory.\n" );

    return EXIT_SUCCESS;

}

int save_field_toHDF5(  gridConfiguration *gridCfg, 
                        saveData *saveDCfg, int t_int,
                        double EB_WAVE[NX][NY][NZ] ){

    if( t_int % (t_save * (int)PERIOD) == 0 && t_int != 0){
        
        size_t ii, jj, kk;
        /*Char values as directions to the correct folder*/
        char filename_hdf5[PATH_MAX], dSet_name[PATH_MAX];
        double (*data2save)[NY/2][NZ/2]     = calloc(NX/2, sizeof *data2save);

        sprintf(filename_hdf5,"%s/%s/%s", projectPath, foldername, file_hdf5);

        // abs(E)
    // prepare array for that
#pragma omp parallel for collapse(3) default(shared) private(ii,jj,kk)
        for (ii=0 ; ii<NX ; ii+=2) {
            for (jj=0 ; jj<NY ; jj+=2) {
                for (kk=0 ; kk<NZ ; kk+=2) {
                    data2save[(ii/2)][(jj/2)][(kk/2)] = 
                        sqrt( pow(EB_WAVE[ii+1][jj  ][kk  ],2) 
                            +pow(EB_WAVE[ii  ][jj+1][kk  ],2) 
                            +pow(EB_WAVE[ii  ][jj  ][kk+1],2) );
                }
            }
        }

        delete_ant2save( gridCfg, data2save );
        
        //Append the name of the files
        sprintf( dSet_name, "E_abs__tint%05d", t_int );
        printf( "status of writeMyHDF_v4: %d\n", writeMyHDF_v4( NX/2, NY/2, NZ/2, filename_hdf5, dSet_name, data2save ) ) ;

        free( data2save );
    }

    return EXIT_SUCCESS;
}

