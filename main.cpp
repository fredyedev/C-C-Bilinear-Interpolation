#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lodepng.cpp"
//#include <iostream>
void freeImage2D(unsigned char **mat);
unsigned char *convert2Dto1D(unsigned char **imgR, unsigned char **imgG,unsigned char **imgB, unsigned int h, unsigned int w);
unsigned char *convert2Dto1D(unsigned char **img2D, unsigned int h, unsigned int w);
unsigned char **getChannel2D(unsigned char *img1D, unsigned int h, unsigned int w,unsigned int index);
unsigned char *readPNG(const char* filename, unsigned int &width, unsigned int &height,unsigned int &bitdepth, unsigned int &bitsXpixel, unsigned int &channels,unsigned int &isGrey, unsigned int &haveAlpha);
int savePNG(const char *filename, unsigned char *img, int width, int height);
double cubicInterpolate (double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
unsigned char **createMatrix(int h,int w) {
	unsigned char **mat;
	int i;
	mat    = (unsigned char **) malloc( h*sizeof(unsigned char *));
    if(mat==NULL) return(NULL);
    mat[0] = (unsigned char *) malloc( h*w*sizeof(unsigned char));
    if(mat[0]==NULL) return(NULL);
    for(i=1; i<h; ++i)
        mat[i] = mat[i-1] + w;
    
	return(mat);
}
unsigned char ** bilinear(int w,int h, int f, unsigned char **imgC){
	unsigned char **img2g = createMatrix(f*h,f*w);
    
	for (int i = 0; i < f*h - 1; i++) {
		
		for (int j = 0; j < f*w - 1; j++) {
			
			int pj = (int)(j / f);
			int pi = (int)(i / f);

			if (pj >= w - 1 || pi >= h - 1) break;

			double fj1 = (double)j / (double)f - (double)pj;
			double fj2 = 1 - fj1;
			double fi1 = (double)i / (double)f - (double)pi;
			double fi2 = 1 - fi1;

			double w1 = fj2*fi2;
			double w2 = fj1*fi2;
			double w3 = fj2*fi1;
			double w4 = fj1*fi1;

			double p1 = (double)imgC[pi][pj];
			double p2 = (double)imgC[pi][pj+1];
			double p3 = (double)imgC[pi+1][pj];
			double p4 = (double)imgC[pi+1][pj+1];
			
			img2g[i][j] = w1*p1 + w2*p2 + w3*p3 + w4*p4;
		}
	}
	return img2g;
}

int main(int argc, char *argv[]) /*list the chunks*/ {
	
    if(argc < 2) {
    	
        printf("De el nombre de la imagen.\n");
        return 0;
        
    }
    
    const char* filename = argv[1];
    char nameout[50];
    int f=atoi(argv[2]);
    

    unsigned char *image1D;
    unsigned int    w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha;

    // Lectura de la imagen
    image1D = readPNG(filename, w, h, bitdepth, bitsXpixel, channels, isGrey, haveAlpha);

    printf("Dimensions:      %d x %d\n", h, w);
    printf("Depth Bits:      %d\n", bitdepth);
    printf("Bits per pixel:  %d\n", bitsXpixel);
    printf("Channels:        %d\n", channels);
    printf("Is Gray Scale?:  %d\n", isGrey);
    printf("Alpha Channel?:  %d\n", haveAlpha);
	sprintf(nameout, "%dx%s_%dx%d.png ",f, "_bl",h,w);
	
    unsigned char **imgR = getChannel2D(image1D, h, w, 0);
    unsigned char **imgG = getChannel2D(image1D, h, w, 1);
    unsigned char **imgB = getChannel2D(image1D, h, w, 2);
    // imagen redimensionada por el factor f
	unsigned char **imgfR;
    unsigned char **imgfG;
    unsigned char **imgfB;
  
  	//para cada canal hacemos interpolacion bilineal

  	imgfR = bilinear(w,h,f, imgR);
	imgfG = bilinear(w,h,f, imgG);
	imgfB = bilinear(w,h,f, imgB);
	
	
    // Se convierte cada canal a un arreglo 1D por separado para poder grabarlos como PNG
    unsigned char *imgR1 = convert2Dto1D(imgfR, f*h, f*w);
    unsigned char *imgG1 = convert2Dto1D(imgfG, f*h, f*w);
    unsigned char *imgB1 = convert2Dto1D(imgfB, f*h, f*w);
	
    // Se combinan los tres canales en un solo arreglo 1D para formar una imagen a color   
   	unsigned char *imgRGB1 = convert2Dto1D(imgfR, imgfG, imgfB, f*h, f*w);
    savePNG(nameout, imgRGB1, f*w, f*h);

    free(image1D);
    freeImage2D(imgR);
    freeImage2D(imgG);
    freeImage2D(imgB);
    freeImage2D(imgfR);
    freeImage2D(imgfG);
    freeImage2D(imgfB);
    free(imgR1);
    free(imgG1);
    free(imgB1);
    free(imgRGB1);
  
    return(0);
}


// Escritura de un arreglo 1D en una imagen PNG
int savePNG(const char *filename, unsigned char *img, int width, int height)
{
    unsigned char *png;
    size_t pngsize;
    int error = lodepng_encode24(&png, &pngsize, img, width, height);
    if(!error)
    {
        lodepng_save_file(png, pngsize, filename);
    }

    if(error)
        printf("\tError %u al grabar la imagen %s: %s\n", error, filename,
                lodepng_error_text(error));

    printf("Termina %s\n", filename);
    free(png);
    return(error);
}


unsigned char *readPNG(const char* filename, unsigned int &width, unsigned int &height,
             unsigned int &bitdepth, unsigned int &bitsXpixel, unsigned int &channels,
             unsigned int &isGrey, unsigned int &haveAlpha) {
  unsigned error;
  unsigned char* image;
  unsigned char* png = 0;
  size_t pngsize;
  LodePNGState state;

  lodepng_state_init(&state);

  error = lodepng_load_file(&png, &pngsize, filename);
  if(!error) error = lodepng_decode(&image, &width, &height, &state, png, pngsize);
  if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
  free(png);

  LodePNGColorMode& color = state.info_png.color;

  bitdepth   = color.bitdepth;
  bitsXpixel = lodepng_get_bpp(&color);
  channels   = lodepng_get_channels(&color);
  isGrey     = lodepng_is_greyscale_type(&color);
  haveAlpha  = lodepng_can_have_alpha(&color);

  lodepng_state_cleanup(&state);
  return(image);
}

// A partir de una imagen img1D codificada como mediante un arreglo 1D, esta funcion
// devuelve un arreglo 2D que tiene la componente de color indicada por index:
// index = 0  para recuperar el canal rojo
// index = 1  para recuperar el canal verde
// index = 2  para recuperar el canal azul
unsigned char **getChannel2D(unsigned char *img1D, unsigned int h, unsigned int w,
                             unsigned int index)
{
    unsigned int    i, j, k, l;
    unsigned int    channels = 3;
    unsigned char **mat;

    // Reservamos memoria
    mat    = (unsigned char **) malloc( h*sizeof(unsigned char *));
    if(mat==NULL) return(NULL);
    mat[0] = (unsigned char *) malloc( h*w*sizeof(unsigned char));
    if(mat[0]==NULL) return(NULL);
    for(i=1; i<h; ++i)
        mat[i] = mat[i-1] + w;

    // Lectura de los datos
    l = (channels + 1);
    for(i=0; i<h; i++) {
        k = i*w*(channels+1);
        for(j=0; j<w; j++) {
            mat[i][j] = img1D[j*l + index + k];
        }
    }

    return(mat);
}

// Convierte un arraglo 2D que tiene la informacion de una imagen a un arreglo 1D
// como se requiere para poderlo grabarlo como imagen PNG en escala de grises
unsigned char *convert2Dto1D(unsigned char **img2D, unsigned int h, unsigned int w)
{
    unsigned char *img1D = (unsigned char *) malloc(sizeof(unsigned char)*h*w*3);
    unsigned int   i, j, k, l;
    unsigned char  val;

    for(i=0; i<h; i++) {
        k = 3*w*i;
        for(j=0; j<w; j++) {
            val = (unsigned char) img2D[i][j];
            l   = 3*j + k;
            img1D[l]   = val;
            img1D[l+1] = val;
            img1D[l+2] = val;
        }
    }
    return(img1D);
}

// Combina la informacion de tres arreglos 2D para formar un arreglo 1D para una imagen a color
unsigned char *convert2Dto1D(unsigned char **imgR, unsigned char **imgG,
                             unsigned char **imgB, unsigned int h, unsigned int w)
{
    unsigned char *img1D = (unsigned char *) malloc(sizeof(unsigned char)*h*w*3);
    unsigned int   i, j, k, l;

    for(i=0; i<h; i++) {
        k = 3*w*i;
        for(j=0; j<w; j++) {
            l   = 3*j + k;
            img1D[l]   = imgR[i][j];
            img1D[l+1] = imgG[i][j];
            img1D[l+2] = imgB[i][j];
        }
    }
    return(img1D);
}


// Libera la memoria del arreglo bidimensional
void freeImage2D(unsigned char **mat) {
    free(mat[0]);
    free(mat);
}

  	
	/*
    for(int i=0; i<2*h; i++) {

        for(int j=0; j<2*w; j++) {
        	
        	img2g[i][j] = (unsigned char) 0;
        	
        }
        
	}*/
/*	
	for(int i=0; i<h; i++) {

        for(int j=0; j<w; j++) {
        	
        	img2g[2*i][2*j] = imgR[i][j];
        	
        }
        
	}
	
	for(int i=0; i<h; i++) {

        for(int j=0; j<w; j++) {
        	
        	if(i!=h-1) {
        		
				img2g[2*i+1][2*j] = (unsigned char) (imgR[i][j]+imgR[i+1][j])/2;
				
        	}
    	}
	}
	
	for(int i=0; i<h; i++) {

        for(int j=0; j<w; j++) {
        	
        	if(j!=w-1) {
        		
				img2g[2*i][2*j+1] = (unsigned char) (imgR[i][j]+imgR[i][j+1])/2;
				
        	}
    	}
	}
	
	for(int i=0; i<h; i++) {

        for(int j=0; j<w; j++) {
        	
        	if(j!=w-1) {
        		
				img2g[2*i+1][2*j+1] = (unsigned char) (imgR[i][j]+imgR[i][j+1])/2;
				
        	}
    	}
	}
	*/
	/*
	for(int i=0; i<h; i++) {

        for(int j=0; j<w; j++) {
        	
        	img2g[i][2*j] = imgR[i][j];
        	
        }
        
	}*/

