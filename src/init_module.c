#include "init_module.h"


void control_init(  gridConfiguration *gridCfg, 
                    beamAntennaConfiguration *beamCfg,
                    saveData *saveDCfg,
                    antennaDetector *antDetect ){
//{{{
        
    /*Initialize System*/
    grid_init( gridCfg, beamCfg, saveDCfg, antDetect );

}//}}}


void grid_init( gridConfiguration *gridCfg, 
                beamAntennaConfiguration *beamCfg,
                saveData *saveDCfg,
                antennaDetector *antDetect ){
    //{{{

    write_JSON_toGrid( gridCfg, beamCfg, saveDCfg, antDetect );

    //Checks that maximum density value is respected
    // if density is larger than this value, FDTD code becomes unstable
    if( ne_0 > PERIOD * 2./5.){
        printf("Density value is too large for code stability. \n");
        printf("Maximum density: %.3f. \n", PERIOD * 2./5.);
        exit(-1);
    }

    //Grid configuration variables computation
    //Set D_ABSORB for the user choose boundary
    if (BOUNDARY == 1){     
        D_ABSORB = (int)(3*PERIOD);
    }else if (BOUNDARY == 2){
        D_ABSORB = 8;
    }

    //Define Z size for reference grid
    if(exc_signal != 6) {
        NZ_REF  = 2*D_ABSORB + (int)PERIOD;
    } else if( exc_signal == 6 ){
        NZ_REF  = 2*D_ABSORB + (int)ant_lenght;
        //NZ_REF  = NZ;
    }

    // dt/dx = 0.5 is commenly used in 2D FDTD codes
    // Note that period refers to the wavelength in the numerical grid and not
    // in the "physical" grid (where one grid cell is equal to one Yee cell).
    // This means that in the physical grid, the wavelength is period/2, thus
    // in the equations we have to use period/2 for the wavelength.
    DX  = 1./(PERIOD/2);
    DT  = 1./(2.*(PERIOD/2));

    // positions have to be even numbers, to ensure fields are accessed correctly
    if ((ANT_X % 2) != 0)  ++ANT_X;
    if ((ANT_Y % 2) != 0)  ++ANT_Y;
    if ((ANT_Z % 2) != 0)  ++ANT_Z;

    //Timetraces number of columns
    col_for_timetraces = 8;

}//}}}


/*Functions in charge of JSON reading*/
void write_JSON_toGrid( gridConfiguration *gridCfg, 
                        beamAntennaConfiguration *beamCfg,
                        saveData *saveDCfg,
                        antennaDetector *antDetect ){
    //{{{

    /*Read JSON and extract data*/
    char *json_file = read_json();
    int scale;

    if(json_file == NULL){
        printf("JSON file doesn't exists.");
        return;
    }

    //Parse the JSON string
    cJSON *json = cJSON_Parse(json_file);
    if(json == NULL){
        printf("Error at parse JSON string");
        free(json_file);
        return;
    }

    //Extract data from JSON and save in struct
    //Save folder info
    cJSON *Main_Project = cJSON_GetObjectItemCaseSensitive(json, "Main_Project");   //Main Project path
    if( cJSON_IsString(Main_Project) && (Main_Project->valuestring != NULL) ){
        projectPath = strdup(Main_Project->valuestring);
    }

    cJSON *FolderName = cJSON_GetObjectItemCaseSensitive(json, "Case_foldername");       //Simulation folder name
    if( cJSON_IsString(FolderName) && (FolderName->valuestring != NULL) ){
        foldername = strdup(FolderName->valuestring);
    }
    
    cJSON *Filename_HDF5 = cJSON_GetObjectItemCaseSensitive(json, "filename_hdf5");   //filename hdf5
    if( cJSON_IsString(Filename_HDF5) && (Filename_HDF5->valuestring != NULL) ){
        file_hdf5 = strdup(Filename_HDF5->valuestring);
    }

    cJSON *Filename_TimeTrace = cJSON_GetObjectItemCaseSensitive(json, "filename_timetraces");   //filename datatraces
    if( cJSON_IsString(Filename_TimeTrace) && (Filename_TimeTrace->valuestring != NULL) ){
        file_trace = strdup(Filename_TimeTrace->valuestring);
    }

    cJSON *t_save_f = cJSON_GetObjectItemCaseSensitive(json, "data_save_frequency");   //time step for data saving
    if( cJSON_IsNumber(t_save_f) ){
        t_save = t_save_f->valueint;
    }

    /*Grid configuration values*/
    cJSON *item_scale = cJSON_GetObjectItemCaseSensitive(json, "scale");   //scale factor
    if( cJSON_IsNumber(item_scale) ){
        scale = item_scale->valueint;
    }

    cJSON *item_period = cJSON_GetObjectItemCaseSensitive(json, "period");   //wave period
    if( cJSON_IsNumber(item_period) ){
        PERIOD = item_period->valuedouble;
        PERIOD = PERIOD * scale;
    }

    cJSON *item_Nx = cJSON_GetObjectItemCaseSensitive(json, "Grid_size_Nx");   //grid size in x
    if( cJSON_IsNumber(item_Nx) ){
        NX = item_Nx->valueint;
        NX = NX * scale;
    }

    cJSON *item_Ny = cJSON_GetObjectItemCaseSensitive(json, "Grid_size_Ny");   //grid size in y 
    if( cJSON_IsNumber(item_Ny) ){
        NY = item_Ny->valueint;
        NY = NY * scale;
    }

    cJSON *item_Nz = cJSON_GetObjectItemCaseSensitive(json, "Grid_size_Nz");   //grid size in z
    if( cJSON_IsNumber(item_Nz) ){
        NZ = item_Nz->valueint;
        NZ = NZ * scale;
    }

    cJSON *item_tend = cJSON_GetObjectItemCaseSensitive(json, "t_end");   //plasma density option
    if( cJSON_IsNumber(item_tend) ){
        T_END = item_tend->valueint;
        T_END = T_END * PERIOD;
    }

    cJSON *item_B0 = cJSON_GetObjectItemCaseSensitive(json, "B0_profile");   //external magnetic field option
    if( cJSON_IsNumber(item_B0) ){
        B0_profile = item_B0->valueint;
    }

    cJSON *item_B0_value = cJSON_GetObjectItemCaseSensitive(json, "B0_magnitude");   //plasma density option
    if( cJSON_IsNumber(item_B0_value) ){
        B0_value = item_B0_value->valuedouble;
    }

    cJSON *item_ne = cJSON_GetObjectItemCaseSensitive(json, "ne_profile");   //plasma density option
    if( cJSON_IsNumber(item_ne) ){
        ne_profile = item_ne->valueint;
    }

    cJSON *item_ne_value = cJSON_GetObjectItemCaseSensitive(json, "ne_0");   //plasma density option
    if( cJSON_IsNumber(item_ne_value) ){
        ne_0 = item_ne_value->valuedouble;
    }

    cJSON *item_boundary = cJSON_GetObjectItemCaseSensitive(json, "Boundary_Method");   //boundary option
    if( cJSON_IsNumber(item_boundary) ){
        BOUNDARY = item_boundary->valueint;
    }

    cJSON *item_dBoundary = cJSON_GetObjectItemCaseSensitive(json, "UPML_size");   //size absorb boundary
    if( cJSON_IsNumber(item_dBoundary) ){
        D_ABSORB = item_dBoundary->valueint;
        D_ABSORB = D_ABSORB + 2;
    }  

    /*Antenna injector configuration values*/
    cJSON *item_ant_x = cJSON_GetObjectItemCaseSensitive(json, "Antenna_Pos_x");   //Antenna x position
    if( cJSON_IsNumber(item_ant_x) ){
        ANT_X = item_ant_x->valuedouble;
    }

    cJSON *item_ant_y = cJSON_GetObjectItemCaseSensitive(json, "Antenna_Pos_y");   //Antenna y position
    if( cJSON_IsNumber(item_ant_y) ){
        ANT_Y = item_ant_y->valuedouble;
    }

    cJSON *item_ant_z = cJSON_GetObjectItemCaseSensitive(json, "Antenna_Pos_z");   //Antenna z position
    if( cJSON_IsNumber(item_ant_z) ){
        ANT_Z = item_ant_z->valuedouble;
    }

    cJSON *item_antAngleZX = cJSON_GetObjectItemCaseSensitive(json, "Antenna_Angle_zx");   //zx plane angle antenna
    if( cJSON_IsNumber(item_antAngleZX) ){
        antAngle_zx = item_antAngleZX->valueint;
    }

    cJSON *item_antAngleZY = cJSON_GetObjectItemCaseSensitive(json, "Antenna_Angle_zy");   //zy plane angle antenna
    if( cJSON_IsNumber(item_antAngleZY) ){
        antAngle_zy = item_antAngleZY->valueint;
    }

    cJSON *item_exc_signal = cJSON_GetObjectItemCaseSensitive(json, "Wave_source_type");   //Signal source option
    if( cJSON_IsNumber(item_exc_signal) ){
        exc_signal = item_exc_signal->valueint;
    }

    cJSON *item_rampUpMethod = cJSON_GetObjectItemCaseSensitive(json, "RampUp_Method");   //RampUp method option
    if( cJSON_IsNumber(item_rampUpMethod) ){
        rampUpMethod = item_rampUpMethod->valueint;
    }

    cJSON *item_ant_w0x = cJSON_GetObjectItemCaseSensitive(json, "Antena_w0x");   //scale factor
    if( cJSON_IsNumber(item_ant_w0x) ){
        ant_w0x = item_ant_w0x->valuedouble;
    }

    cJSON *item_ant_w0y = cJSON_GetObjectItemCaseSensitive(json, "Antena_w0y");   //scale factor
    if( cJSON_IsNumber(item_ant_w0y) ){
        ant_w0y = item_ant_w0y->valuedouble;
    }

    cJSON *item_z2waist = cJSON_GetObjectItemCaseSensitive(json, "z2waist");   //
    if( cJSON_IsNumber(item_z2waist) ){
        z2waist = item_z2waist->valuedouble;
        z2waist = z2waist * .0;            // .2/l_0*period = -298.87
    }

    /*Antenna Detector Input values*/
    cJSON *item_antDetect = cJSON_GetObjectItemCaseSensitive(json, "Detector_Antenna");   //Activate Antenna
    if( cJSON_IsNumber(item_antDetect) ){
        activate_antDetect1D = item_antDetect->valueint;
    }

    /*Helical Antenna*/
    cJSON *item_Antenna_type = cJSON_GetObjectItemCaseSensitive(json, "Antenna_type");   //Helical model
    if( cJSON_IsNumber(item_Antenna_type) ){
        ant_type = item_Antenna_type->valueint;
    }

    cJSON *item_Antenna_lenght = cJSON_GetObjectItemCaseSensitive(json, "Antenna_lenght");   //Helical model
    if( cJSON_IsNumber(item_Antenna_lenght) ){
        ant_lenght = item_Antenna_lenght->valueint;
    }

    cJSON *item_Radius = cJSON_GetObjectItemCaseSensitive(json, "Radius");   //
    if( cJSON_IsNumber(item_Radius) ){
        ant_radius = item_Radius->valuedouble;
    }

    cJSON *item_J_amp = cJSON_GetObjectItemCaseSensitive(json, "Currrent_Amp");   //
    if( cJSON_IsNumber(item_J_amp) ){
        J_amp = item_J_amp->valuedouble;
    }

    //clean up
    cJSON_Delete(json);
    free(json_file);

}//}}}


char *read_json(){
    //{{{

    FILE *file = fopen("input_FOCAL.json", "rb");
    if (file == NULL) {
        perror("Error openning file.");
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);

    char *json_data = (char *)malloc(length + 1);
    if (json_data == NULL) {
        perror("Error at allocating memory");
        fclose(file);
        return NULL;
    }

    fread(json_data, 1, length, file);
    json_data[length] = '\0';

    fclose(file);
    return json_data;

}//}}}


/*Configuration print on console*/
void print_systemConfiguration(gridConfiguration *gridCfg, beamAntennaConfiguration *beamCfg ){
    //{{{

    // print some info to console
    printf("----------System Configuration Parameters----------\n");
    printf( "Nx = %d, Ny = %d, Nz = %d\n", NX, NY, NZ );
    printf( "period = %d\n", (int)(PERIOD) );
    printf( "d_absorb = %d\n", D_ABSORB );
    printf( "t_end = %d\n", (int)(T_END) );
    printf( "antAngle_zx = %.2f, antAngle_zy = %.2f\n", antAngle_zx, antAngle_zy );
    printf( "ant_w0x = %.2f, ant_w0y = %.2f\n", ant_w0x, ant_w0y ); 
    printf( "ant_x = %d, ant_y = %d, ant_z = %d\n", ANT_X, ANT_Y, ANT_Z );
    printf( "Boundary condition set to '%d'\n", BOUNDARY );
    printf( "Courant number = %.2f. \n", DT/DX);

}//}}}

