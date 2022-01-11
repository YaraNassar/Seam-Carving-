#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void calc_energy(struct rgb_img *im, struct rgb_img **grad){ 
    int h = im->height;
    int w = im->width;
    *grad = (struct rgb_img *)malloc(sizeof(struct rgb_img));
    (*grad)->height = h;
    (*grad)->width = w;
    (*grad)->raster = (uint8_t *)malloc(3 * h * w);

    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            int left = j - 1;
            int right = j + 1;
            int top = i - 1;
            int bot = i + 1;
            if(left == -1){
                left = w - 1;
            }
            if(right == w){
                right = 0;
            }
            if(top == -1){
                top = h - 1;
            }
            if(bot == h){
                bot = 0;
            }

            double R_r = get_pixel(im, i, right, 0);
            double G_r = get_pixel(im, i, right, 1);
            double B_r = get_pixel(im, i, right, 2);

            double R_l = get_pixel(im, i, left, 0);
            double G_l = get_pixel(im, i, left, 1);
            double B_l = get_pixel(im, i, left, 2);

            double R_x = R_r - R_l;
            double G_x = G_r - G_l;
            double B_x = B_r - B_l;

            double delta_x_sq = pow(R_x, (double)2) + pow(G_x, (double)2) + pow(B_x, (double)2);

            double R_t = get_pixel(im, top, j, 0);
            double G_t = get_pixel(im, top, j, 1);
            double B_t = get_pixel(im, top, j, 2);

            double R_b = get_pixel(im, bot, j, 0);
            double G_b = get_pixel(im, bot, j, 1);
            double B_b = get_pixel(im, bot, j, 2);

            double R_y = R_b - R_t;
            double G_y = G_b - G_t;
            double B_y = B_b - B_t;

            double delta_y_sq = pow(R_y, (double)2) + pow(G_y, (double)2) + pow(B_y, (double)2);

            double energy = sqrt(delta_x_sq + delta_y_sq);
            uint8_t e_p = (uint8_t)(energy / (double)10);
            set_pixel((*grad), i, j, e_p, e_p, e_p);
        }
    }
}
double min1(double x, double y, double z){
    double m;
    if(x < y){
        m = x;
    }else{
        m = y;
    }
    if(z < m){
        m = z;
    }
    return m;
}

double min2(double x, double y){
    double m;
    if(x < y){
        m = x;
    }else{
        m = y;
    }
    return m;
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){
    int height = grad->height;
    int width = grad->width;
    (*best_arr) = (double *)malloc(height * width * sizeof(double));
    double *arr = (double *)malloc(height * width * sizeof(double));
    for(int j = 0; j < height * width; j++){
        double num = (double)grad->raster[3*j];
        arr[j] = num;
    }
    double m;
    int i;
    for(i = 0; i < width; i++){
        (*best_arr)[i] = arr[i];
    }
    while(i < height * width){
        if(i % width == 0){
            m = min2((*best_arr)[i - width], (*best_arr)[i - width + 1]);  
        }else if((i + 1) % width == 0){
            m = min2((*best_arr)[i - width], (*best_arr)[i - width - 1]);
        }else{
            m = min1((*best_arr)[i - width], (*best_arr)[i - width + 1], (*best_arr)[i - width - 1]);
        } 
        (*best_arr)[i] = arr[i] + m;
        i++;
    }
    free(arr);

}


int min_list(double *best, int width, int height){
    int i = (height * width) - width;
    double min = best[i];
    for(int j = (height * width) - width + 1; j < height * width; j++){
        if(best[j] < min){
            min = best[j];
            i = j;
        }
    }
    return i;
}


void recover_path(double *best, int height, int width, int **path){
    (*path) = (int *)malloc(height * sizeof(int));
    int i = min_list(best, width, height);
    
    int level = height - 1;
    (*path)[level] = i % width;
    for(level = height - 2; level >= 0; level = level - 1){
        if(i % width == 0){
            if(best[i - width] < best[i - width + 1]){
                i = i - width;
            }else{
                i = i - width + 1;
            }
        }else if((i + 1) % width == 0){
            if(best[i - width] < best[i - width - 1]){
                i = i - width;
            }else{
                i = i - width - 1;
            }
        }else{
            int m;
            if(best[i - width + 1] < best[i - width]){
                m = i - width + 1 ;
            }else{
                m = i - width;
            }
            if(best[i - width - 1] < best[m]){
                m = i - width -1;
            }
            i = m;
        }
        (*path)[level] = i % width;
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path){
    int height = src->height;
    int width = src->width;
    create_img(dest, height, width - 1);
    int i = 0;
    for(int level = 0; level < height; level++){
        for(; i < 3 * ((level*width) + path[level]); i++){
            (*dest)->raster[i - 3*level] = src->raster[i];
        }
        i++;
        i++;
        i++;
        for(; i < 3*(width * (1 + level)); i++){
            (*dest)->raster[i - 3*(level + 1)] = src->raster[i];
        }
    }

}
 