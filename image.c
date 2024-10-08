#include "lodepng.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



char* load_png_file  (const char *filename, int *width, int *height) {
	unsigned char *image = NULL;
	int error = lodepng_decode32_file(&image, width, height, filename);
	if (error) {
		printf("error %u: %s\n", error, lodepng_error_text(error));
		return NULL;
	}

	return (image);
}

void Prewitt_Sobel(unsigned char* image, int width, int height) {
    int i, j, t;
    unsigned char* edge_image=malloc(width*height*4*sizeof(unsigned char));
    for (int i=1;i<height-1;i++){
        for (int j=1;j<width-1;j++){
            double gradient_x = (double) (((double) image[4*((i+1)*width+(j-1))]+1.5*image[4*((i+1)*width+j)]+image[4*((i+1)*width+(j+1))])-
                              ((double) image[4*((i-1)*width+(j-1))]+1.5*image[4*((i-1)*width+j)]+image[4*((i-1)*width+(j+1))]));

            double gradient_y = (double) (((double) image[4*((i-1)*width+(j+1))]+1.5*image[4*(i*width+(j+1))]+image[4*((i+1)*width+(j+1))])-
                              ((double) image[4*((i-1)*width+(j-1))]+1.5*image[4*(i*width+(j-1))]+image[4*((i+1)*width+(j-1))]));

            int gradient = (int)sqrt((double)(gradient_x*gradient_x+gradient_y*gradient_y));
            int k;
            for (k=0; k<3; k++){
                edge_image[4*(i*width+j)+k]=(gradient>15)?255:image[4*(i*width+j)+k];
            }
            edge_image[4*(i*width+j)+3]=image[4*(i*width+j)+3];
        }
    }
    for (int i=0;i<height*width*4;i++){
        image[i]=edge_image[i];
    }
    free(edge_image);
}

typedef struct Set{
    int representative;
    int rank;
}Set;

void Make_Set(Set* sets, int v){
    sets[v].representative=v;
    sets[v].rank=0;
    return;
}

int Find_Set(Set* sets, int v){
    if(v!=sets[v].representative){
        sets[v].representative=Find_Set(sets,sets[v].representative);
    }
    return sets[v].representative;
}

void Union(Set* sets, int i, int j){
    int root1=Find_Set(sets, i);
    int root2=Find_Set(sets, j);
    if(root1!=root2){
        if (sets[root1].rank<sets[root2].rank){
            sets[root1].representative=root2;
        }else if(sets[root1].rank>sets[root2].rank){
            sets[root2].representative=root1;
        }else{
            sets[root2].representative=root1;
            sets[root1].rank++;
        }
    }
    return;
}



int main(){
    unsigned ws, hs, wh, hh;
    int i, j;
    const char* filename1="skull.png";
    const char* filename2="hand.png";
    unsigned char* image_skull=NULL;
    unsigned char* image_hand=NULL;
    unsigned error1=lodepng_decode32_file(&image_skull, &ws, &hs, filename1);
    unsigned error2=lodepng_decode32_file(&image_hand, &wh, &hh, filename2);
    unsigned char* grey_image_skull=malloc(4*ws*hs*sizeof(unsigned char));
    unsigned char* grey_image_hand=malloc(4*wh*hh*sizeof(unsigned char));
    for (unsigned i=0;i<4*ws*hs;i+=4){
        grey_image_skull[i]=(image_skull[i]+image_skull[i+1]+image_skull[i+2])/3;
        grey_image_skull[i+1]=(image_skull[i]+image_skull[i+1]+image_skull[i+2])/3;
        grey_image_skull[i+2]=(image_skull[i]+image_skull[i+1]+image_skull[i+2])/3;
        grey_image_skull[i+3]=255;
    }

    for (unsigned i=0;i<4*wh*hh;i+=4){
        grey_image_hand[i]=(image_hand[i]+image_hand[i+1]+image_hand[i+2])/3;
        grey_image_hand[i+1]=(image_hand[i]+image_hand[i+1]+image_hand[i+2])/3;
        grey_image_hand[i+2]=(image_hand[i]+image_hand[i+1]+image_hand[i+2])/3;
        grey_image_hand[i+3]=255;
    }

    unsigned char* blur_image_skull = (unsigned char*) malloc(4*ws*hs*sizeof(unsigned char));
    unsigned char* blur_image_hand = (unsigned char*) malloc(4*wh*hh*sizeof(unsigned char));

    int radius=3;
    float sigma=6.0;
    int core_size=2*radius+1;
    float core[core_size*core_size];
    float core_sum=0;
    for (int i=-radius;i<=radius;i++){
        for (int j=-radius;j<=radius;j++){
            int idx = (i+radius)*core_size+(j+radius);
            core[idx]=(float) expf(-(i*i+j*j)/(2.0*sigma*sigma))/(2.0*M_PI*sigma*sigma);
            core_sum+=core[idx];
        }
    }
    for (i=0;i<core_size*core_size;i++){
        core[i]=core[i]/core_sum;
    }
    int x, y;
    for (y=0; y<hs; y++){
        for (x=0; x<ws; x++){
            float grey=0;
            for (int i = -radius; i <= radius; i++){
                for (int j=-radius; j<=radius; j++){
                    int nx=x+j;
                    int ny=y+i;
                    if (nx<0){
                            nx=0;
                    }
                    if (nx>=ws){
                            nx=ws-1;
                    }
                    if (ny<0){
                            ny=0;
                    }
                    if (ny>=hs){
                            ny=hs-1;
                    }

                    int offset=(ny*ws+nx)*4;
                    float core_value=core[(i+radius)*core_size+(j+radius)];
                    grey+=image_skull[offset]*core_value;
                }
            }

            int current_offset=(y*ws+x)*4;
            blur_image_skull[current_offset]=255;
            blur_image_skull[current_offset+1]=255;
            blur_image_skull[current_offset+2]=255;
            blur_image_skull[current_offset+3]=255;
            if (grey<255){
                blur_image_skull[current_offset] = (unsigned char) grey;
                blur_image_skull[current_offset+1] = (unsigned char) grey;
                blur_image_skull[current_offset+2] = (unsigned char) grey;
                blur_image_skull[current_offset+3] = (unsigned char) grey;
            }
        }
    }

    core_sum=0;
    for (int i=-radius;i<=radius;i++){
        for (int j=-radius;j<=radius;j++){
            int idx = (i+radius)*core_size+(j+radius);
            core[idx]=(float) expf(-(i*i+j*j)/(2.0*sigma*sigma))/(2.0*M_PI*sigma*sigma);
            core_sum+=core[idx];
        }
    }
    for (i=0;i<core_size*core_size;i++){
        core[i]=core[i]/core_sum;
    }
    for (y=0; y<hh; y++){
        for (x=0; x<wh; x++){
            float grey=0;
            for (int i = -radius; i <= radius; i++){
                for (int j=-radius; j<=radius; j++){
                    int nx=x+j;
                    int ny=y+i;
                    if (nx<0){
                            nx=0;
                    }
                    if (nx>=ws){
                            nx=wh-1;
                    }
                    if (ny<0){
                            ny=0;
                    }
                    if (ny>=hs){
                            ny=hh-1;
                    }

                    int offset=(ny*wh+nx)*4;
                    float core_value=core[(i+radius)*core_size+(j+radius)];
                    grey+=image_hand[offset]*core_value;
                }
            }

            int current_offset=(y*wh+x)*4;
            blur_image_hand[current_offset]=255;
            blur_image_hand[current_offset+1]=255;
            blur_image_hand[current_offset+2]=255;
            blur_image_hand[current_offset+3]=255;
            if (grey<255){
                blur_image_hand[current_offset] = (unsigned char) grey;
                blur_image_hand[current_offset+1] = (unsigned char) grey;
                blur_image_hand[current_offset+2] = (unsigned char) grey;
                blur_image_hand[current_offset+3] = (unsigned char) grey;
            }
        }
    }

    Prewitt_Sobel(blur_image_skull, ws, hs);
    Prewitt_Sobel(blur_image_hand, wh, hh);

    unsigned error3=lodepng_encode32_file("preskull.png", blur_image_skull, ws, hs);
    unsigned error4=lodepng_encode32_file("prehand.png", blur_image_hand, wh, hh);
    int Vs=hs*ws;
    int Vh=hh*wh;
    Set* setss = (Set*) malloc(Vs*sizeof(Set));
    Set* setsh = (Set*) malloc(Vs*sizeof(Set));
    for (int i=0;i<Vs;i++) Make_Set(setss,i);
    for (int i=0;i<Vh;i++) Make_Set(setsh,i);
    unsigned char* skull = (unsigned char*) malloc(4*ws*hs*sizeof(unsigned char));
    unsigned char* hand = (unsigned char*) malloc(4*wh*hh*sizeof(unsigned char));

    unsigned char delta = 2;

    skull[0]=0;
    skull[1]=0;
    skull[2]=0;
    skull[3]=255;
    for (i=0;i<hs;i++){
        for (j=0;j<ws;j++){
            int pixel0=4*(i*ws+j);
            int vertex0=i*ws+j;
            int pixel1=4*(i*ws+j+1);
            int vertex1=i*ws+j+1;
            int pixel2=4*((i+1)*ws+j);
            int vertex2=(i+1)*ws+j;
            int pixel3=4*((i+1)*ws+j+1);
            int vertex3=(i+1)*ws+j+1;
            if(pixel0+3>=4*hs*ws || pixel1+3>=4*hs*ws || pixel2+3>=4*hs*ws || pixel3+3>=4*hs*ws){

            } else {
                unsigned char x0=blur_image_skull[pixel0];
                unsigned char x1=-1;
                unsigned char x2=-1;
                unsigned char x3=-1;

                if (j+1<ws){
                    x1 = blur_image_skull[pixel1];
                    if (i+1<hs){
                        x2 = blur_image_skull[pixel2];
                        x3 = blur_image_skull[pixel3];
                    }
                } else{
                    if (i+1<hs){
                        x2 = blur_image_skull[pixel2];
                    }
                }


                if((x0-x1)*(x0-x1)<=delta && x1!=-1){
                    Union(setss,vertex0,vertex1);
                }
                if((x0-x2)*(x0-x2)<=delta && x2!=-1){
                    Union(setss,vertex0,vertex2);
                }
                if((x0-x3)*(x0-x3)<=delta && x3!=-1){
                    Union(setss,vertex0,vertex3);
                }
            }
        }
    }

    hand[0]=0;
    hand[1]=0;
    hand[2]=0;
    hand[3]=255;
    for (i=0;i<hh;i++){
        for (j=0;j<wh;j++){
            int pixel0=4*(i*wh+j);
            int vertex0=i*wh+j;
            int pixel1=4*(i*wh+j+1);
            int vertex1=i*wh+j+1;
            int pixel2=4*((i+1)*wh+j);
            int vertex2=(i+1)*wh+j;
            int pixel3=4*((i+1)*wh+j+1);
            int vertex3=(i+1)*wh+j+1;
            if(pixel0+3>=4*hh*wh || pixel1+3>=4*hh*wh || pixel2+3>=4*hh*wh || pixel3+3>=4*hh*wh){
                //continue;
            } else {
                unsigned char x0=blur_image_hand[pixel0];
                unsigned char x1=-1;
                unsigned char x2=-1;
                unsigned char x3=-1;

                if (j+1<wh){
                    x1 = blur_image_hand[pixel1];
                    if (i+1<hh){
                        x2 = blur_image_hand[pixel2];
                        x3 = blur_image_hand[pixel3];
                    }
                } else{
                    if (i+1<hh){
                        x2 = blur_image_hand[pixel2];
                    }
                }


                if((x0-x1)*(x0-x1)<=delta && x1!=-1){
                    Union(setsh,vertex0,vertex1);
                }
                if((x0-x2)*(x0-x2)<=delta && x2!=-1){
                    Union(setsh,vertex0,vertex2);
                }
                if((x0-x3)*(x0-x3)<=delta && x3!=-1){
                    Union(setsh,vertex0,vertex3);
                }
            }
        }
    }


    int R=0, G=0, B=0;
    int* roots=malloc(Vs*sizeof(int));
    int** components=malloc(Vs*sizeof(int*));
    for (i=0;i<Vs;i++){
            components[i]= (int**) malloc(sizeof(int));
    }
    for (i=0;i<Vs;i++){
        int root=Find_Set(setss,i);
        roots[i]=root;
        components[root] = realloc(components[root], ((components[root][0])+2)*sizeof(int));
        components[root][0]++;
        components[root][components[root][0]]=i;
    }
    for (i=0;i<Vs;i++){
        if(roots[i]==i){
            R=(R+177)%256;
            G=(G+215)%256;
            B=(B+123)%256;
            for (j=1;j<=components[i][0];j++){
                int current=components[i][j]+1;
                if ((4*current+3)>=4*Vs){
                } else {
                    skull[4*current]=R;
                    skull[4*current+1]=G;
                    skull[4*current+2]=B;
                    skull[4*current+3]=255;
                }
            }
        }
    }
    for(i=0;i<Vs;i++){
            free(components[i]);
    }
    free(components);
    free(roots);



    R=0;
    G=0;
    B=0;
    roots=malloc(Vh*sizeof(int));
    components=malloc(Vh*sizeof(int*));
    for (i=0;i<Vh;i++){
            components[i] = (int**) malloc(sizeof(int));
    }
    for (i=0;i<Vh;i++){
        int root=Find_Set(setsh,i);
        roots[i]=root;
        components[root] = realloc(components[root], ((components[root][0])+2)*sizeof(int));
        components[root][0]++;
        components[root][components[root][0]]=i;
    }
    for (int i=0;i<Vh;i++){
        if(roots[i]==i){
            R=(R+177)%256;
            G=(G+215)%256;
            B=(B+123)%256;
            for (j=1;j<=components[i][0];j++){
                int current=components[i][j]+1;
                if ((4*current+3)>=4*Vh){
                } else {
                    hand[4*current]=R;
                    hand[4*current+1]=G;
                    hand[4*current+2]=B;
                    hand[4*current+3]=255;
                }
            }
        }
    }
    for(i=0;i<Vh;i++){
            free(components[i]);
    }
    free(components);
    free(roots);


    unsigned error5=lodepng_encode32_file("outskull.png", skull, ws, hs);
    unsigned error6=lodepng_encode32_file("outhand.png", hand, wh, hh);

    free(setss);
    free(setsh);
    free(image_skull);
    free(image_hand);
    free(grey_image_skull);
    free(grey_image_hand);
    free(blur_image_skull);
    free(blur_image_hand);
    free(skull);
    free(hand);
    return 0;
}
