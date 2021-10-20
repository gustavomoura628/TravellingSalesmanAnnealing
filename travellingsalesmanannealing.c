//INCLUDE SDL
#include <SDL2/SDL.h> 
//#include <SDL2/SDL_image.h> 
#include <SDL2/SDL_mouse.h>
#include <SDL2/SDL_rect.h>
#include <SDL2/SDL_render.h>
#include <SDL2/SDL_scancode.h>
#include <SDL2/SDL_timer.h> 
#include <SDL2/SDL_ttf.h>

//INCLUDE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <complex.h>


//SDL Related Functions
void CSDL_RenderDrawRectangle(SDL_Renderer* renderer, int x, int y, int w, int h)
{
	int i;
    SDL_Rect r;
    r.x = x;
    r.y = y;
    r.w = w;
    r.h = h;
    SDL_RenderFillRect(renderer,&r);
}
void CSDL_RenderDrawCircle(SDL_Renderer* renderer, int x, int y, int r)
{
	int i,j;
	j=0;
	i=r;
	int k;
	int f=0;
	for(j=0;i>=j;j++){
		f=j*j+i*i-r*r;
		if(f>0)i--;
		SDL_RenderDrawLine(renderer, x-j, y+i, x+j, y+i);
		SDL_RenderDrawLine(renderer, x-j, y-i, x+j, y-i);
		SDL_RenderDrawLine(renderer, x-i, y+j, x+i, y+j);
		SDL_RenderDrawLine(renderer, x-i, y-j, x+i, y-j);
	}
}
void CSDL_RenderDrawLine(SDL_Renderer* renderer, int x0, int y0, int x1, int y1, int r)
{
    int dx = abs(x0-x1);
    int dy = abs(y0-y1);

	int i,j;
	for(i=-r/2;i<r/2;i++){
        if(dx<dy){
            SDL_RenderDrawLine(renderer, x0+i, y0, x1+i, y1);
        }else{
            SDL_RenderDrawLine(renderer, x0, y0+i, x1, y1+i);
        }
	}
}
void CSDL_RenderDrawRoundLine(SDL_Renderer* renderer, int x0, int y0, int x1, int y1, int r)
{
    CSDL_RenderDrawLine(renderer,x0,y0,x1,y1,r);
    CSDL_RenderDrawCircle(renderer,x0,y0,r/2);
    CSDL_RenderDrawCircle(renderer,x1,y1,r/2);
}
void CSDL_RenderDrawText(SDL_Renderer * renderer, char * textInput, int x, int y, int textSize, int outline_size)
{
    TTF_Font *Buran = TTF_OpenFont("fonts/Buran.ttf", textSize); //this opens a font style and sets a size
    TTF_Font *font_outline = TTF_OpenFont("fonts/Buran.ttf", textSize); 
    TTF_SetFontOutline(font_outline, outline_size); 
    SDL_Color textColor;
    SDL_GetRenderDrawColor(renderer, &textColor.r, &textColor.g, &textColor.b, &textColor.a);
    //= {r, g, b, 0};
    SDL_Color inverseTextColor = {255 - textColor.r, 255 - textColor.g, 255 - textColor.b, 0};
    SDL_Surface* surfaceMessage = TTF_RenderText_Blended(Buran, textInput, textColor); // as TTF_RenderText_Solid could only be used on SDL_Surface then you have to create the surface first
    SDL_Surface* surfaceOutline = TTF_RenderText_Blended(font_outline, textInput, inverseTextColor); // as TTF_RenderText_Solid could only be used on SDL_Surface then you have to create the surface first
    TTF_CloseFont(Buran);
    TTF_CloseFont(font_outline);
    SDL_Rect rect = { outline_size, outline_size, surfaceMessage->w, surfaceMessage->h}; 
    SDL_SetSurfaceBlendMode(surfaceMessage, SDL_BLENDMODE_BLEND);
    SDL_BlitSurface(surfaceMessage, NULL, surfaceOutline, &rect);
    //SDL_Texture* Message = SDL_CreateTextureFromSurface(renderer, surfaceMessage); //now you can convert it into a texture
    SDL_Texture* Outline = SDL_CreateTextureFromSurface(renderer, surfaceOutline); //now you can convert it into a texture
    SDL_FreeSurface(surfaceMessage); 
    SDL_FreeSurface(surfaceOutline); 
    SDL_Rect Outline_rect; //create a rect
    SDL_QueryTexture(Outline, NULL, NULL, &Outline_rect.w, &Outline_rect.h); 
    Outline_rect.x = x-(Outline_rect.w/2);  //controls the rect's x coordinate 
    Outline_rect.y = y-(Outline_rect.h/2); // controls the rect's y coordinte
    SDL_RenderCopy(renderer, Outline, NULL, &Outline_rect); 
    SDL_DestroyTexture(Outline); 	
}

float myrand(float a, float b)
{
    float r = rand()/(float)RAND_MAX;
    return r*fabs(a-b)+fmin(a,b);
}

//HamiltonianCycle data structure that will be used for the annealing process
typedef struct hamiltoniancycle{
    int * order;
    double ** weights;
    int n;
}HC;

void HCDestroy(HC * hc){
    free(hc->order);
}
HC * HCCopy(HC * hc){
    HC * hcc = (HC*)malloc(sizeof(HC));
    hcc->n = hc->n;
    hcc->order = (int*)malloc(hcc->n * sizeof(int));
    hcc->weights = hc->weights;
    memcpy(hcc->order, hc->order, hc->n * sizeof(int));
    return hcc;
}
HC * HCGetSuccessor(HC * hc){
    HC * t = HCCopy(hc);
    int i0 = rand()%(t->n);
    int i1=i0;
    while(i1==i0)i1 = rand()%(t->n);
    int temp = t->order[i0];
    t->order[i0] = t->order[i1];
    t->order[i1] = temp;
    return t;
}
float HCScore(HC * hc){
    float score = 0;
    int i;
    for(i=0;i<hc->n;i++){
        score-=hc->weights[hc->order[i]][hc->order[(i+1)%hc->n]];
    }
    return score;
}

//Simulated annealing functions
float simulatedAnnealingScoreCallback(void * state){
    return HCScore(state);
}
void simulatedAnnealingDestroyCallback(void * state){
    HCDestroy(state);
    free(state);
}
void * simulatedAnnealingSuccessorCallback(void * state){
    return HCGetSuccessor(state);

}
void * simulatedAnnealing(void * state, float temp, void * getSuccessor(void * state), float getScore(void * state), void destroyState(void * state)){
    void * s = state;
    void * t = getSuccessor(s);
    temp = temp < 0 ? 0 : temp;
    float sscore = getScore(s);
    float tscore = getScore(t);
    if(tscore >= sscore){
        destroyState(s);
        s = t;
    }else{
        float loss = fabs(sscore - tscore);
        float probability = exp(-(loss/temp));
        if(rand()/(float)RAND_MAX <= probability){
            destroyState(s);
            s=t;
        }else{
            destroyState(t);
        }

    }
    return s;
}

int main(int argc, char ** argv) 
{ 
    int i,j,k;
    //Initialize srand
    struct timespec t;
    timespec_get(&t,TIME_UTC);
	srand(t.tv_nsec ^ t.tv_sec);

    //Initialize SDL
	if (SDL_Init(SDL_INIT_EVERYTHING) != 0) printf("error initializing SDL: %s\n", SDL_GetError()); 
	 
    //Gets display 0's resolution
    SDL_DisplayMode SDL_SCREENMODE;
    SDL_GetDisplayMode(0,0,&SDL_SCREENMODE);
    const int SCREEN_WIDTH  = SDL_SCREENMODE.w;
    const int SCREEN_HEIGHT = SDL_SCREENMODE.h;
    float FPS = 60;


	SDL_Window* window = SDL_CreateWindow("GAME", // creates a window 
									SDL_WINDOWPOS_CENTERED, 
									SDL_WINDOWPOS_CENTERED, 
									SCREEN_WIDTH, SCREEN_HEIGHT, 0); 
    SDL_SetWindowFullscreen(window,SDL_WINDOW_FULLSCREEN);

	Uint32 render_flags = SDL_RENDERER_ACCELERATED; 
	SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, render_flags); 
	TTF_Init();
    double SCREEN_K = SCREEN_HEIGHT+SCREEN_WIDTH;

    //controls animation loop 
	int close = 0; 
	int pause = 1;

    //Declares input handling variables
	int keyboard[1000];
	int keyboardprev[1000];
	int keyboardonce[1000];
	int keyboarddown[1000];
	int keyboardup[1000];
	int mouse[10];
	int mouseprev[10];
	int mouseonce[10];
	int mousedown[10];
	int mouseup[10];
	unsigned int counter;
	int lastTime = SDL_GetTicks();
	memset(keyboard,0,sizeof keyboard);
	
	int mouseX;
	int mouseY;

    //Declares Coordinate related variables
    double kx = .4758*SCREEN_WIDTH;
    double ky = -.8463*SCREEN_HEIGHT;
    double bx = 0.1;
    double by = SCREEN_HEIGHT/ky+0.1;

    //moving coordinate system:
    double coordmouseclickdownx = 0;
    double coordmouseclickdowny = 0;
    double coordmouseoriginalbx = 0;
    double coordmouseoriginalby = 0;
    int coormousetranslation = 0;

    //Initializes Points
    int npoints = 50;
    double ** points = (double**)malloc(npoints*sizeof(double*));
    for(i=0;i<npoints;i++){
        points[i] = (double*)malloc(sizeof(double)*2);
        points[i][0] = myrand(0,1)+0.05;
        points[i][1] = myrand(0,1)+0.05;
    }

    //Initializes Vertices
    int NumberOfVertices = npoints;
    double ** weights = (double **)malloc(NumberOfVertices*sizeof(double*));
    for(i=0;i<NumberOfVertices;i++){
        weights[i] = (double*)malloc(NumberOfVertices*sizeof(double));
    }

    //Initializes Edges
    int NumberOfEdges = NumberOfVertices*NumberOfVertices;
    for(i=0;i<NumberOfVertices;i++){
        for(j=i+1;j<NumberOfVertices;j++){
                double w = pow(points[i][0]-points[j][0],2) + pow(points[i][1]-points[j][1],2);
                weights[i][j] = w;
                weights[j][i] = w;
        }
    }

    //Initializes the HamiltonianCycle structure that the annealing will be operated on
    HC * hc = (HC*)malloc(sizeof(HC));
    hc->weights = weights;
    hc->order = (int*)malloc(NumberOfVertices*sizeof(int));
    hc->n = npoints;
    for(i=0;i<NumberOfVertices;i++)hc->order[i]=i;

    HC * hcmax = HCCopy(hc);
    int uncrossing = 0;
    int canbeuncrossed = 1;
    HC * hcuncross;


    //Initializes annealing variables
    float temp = 5;
    float coolingRate = 0.997;

	//Animation loop 
	while (!close)
    {
        // Events mangement 
		SDL_Event event; 
		while (SDL_PollEvent(&event)) {
			switch (event.type) { 

			case SDL_QUIT: 
				// handling of close button 
				close = 1; 
				break; 
			case SDL_KEYDOWN:
				keyboard[event.key.keysym.scancode] = 1;
				keyboarddown[event.key.keysym.scancode] = 1;
				break;
			case SDL_KEYUP:
				keyboard[event.key.keysym.scancode] = 0;
				keyboardup[event.key.keysym.scancode] = 1;
				break;
			case SDL_MOUSEBUTTONDOWN:
				mouse[event.button.button]=1;
				mousedown[event.button.button]=1;
				break;
			case SDL_MOUSEBUTTONUP:
				mouse[event.button.button]=0;
				mouseup[event.button.button]=1;
				break;
            case SDL_MOUSEWHEEL:
                if(event.wheel.y > 0) // scroll up
                {
                    //Handles zoom
                    double Sx = (mouseX)/kx-bx;
                    double Sy = (mouseY)/ky-by;
                    kx*=1.1;
                    ky*=1.1;
                    bx = (mouseX)/kx - Sx;
                    by = (mouseY)/ky - Sy;
                }
                else if(event.wheel.y < 0) // scroll down
                {
                    //Handles zoom
                    double Sx = (mouseX)/kx-bx;
                    double Sy = (mouseY)/ky-by;
                    kx*=1/(1.1);
                    ky*=1/(1.1);
                    bx = (mouseX)/kx - Sx;
                    by = (mouseY)/ky - Sy;
                }
                break;
			}
		} 

        //Input management
        //handle keyboardonce variable
		for(i=0;i<1000;i++)
		{
			keyboardonce[i] = 0;
			if(!keyboardprev[i] && keyboard[i])keyboardonce[i] = 1;
			keyboardprev[i] = keyboard[i];
		}
        //handle mouseonce variable
		for(i=0;i<10;i++){
			mouseonce[i] = 0;
			if(!mouseprev[i] && mouse[i])mouseonce[i]=1;
			mouseprev[i]=mouse[i];
		}
        SDL_GetMouseState(&mouseX,&mouseY);

        //handle actions
		if(keyboardonce[SDL_SCANCODE_RETURN])
		{
		}
		if(keyboard[SDL_SCANCODE_R])
		{
            //resets coordinate system
            kx = 100;
            ky = -100;
            bx = SCREEN_WIDTH/2.0/kx;
            by = SCREEN_HEIGHT/2.0/ky;
        }

		if(keyboardonce[SDL_SCANCODE_K])
		{
		}
		if(keyboardonce[SDL_SCANCODE_J])
		{
		}

		if(keyboardonce[SDL_SCANCODE_N])
		{
		}
		if(keyboardonce[SDL_SCANCODE_P])
		{
		}
		if(keyboardonce[SDL_SCANCODE_L])
		{
		}
		if(keyboardonce[SDL_SCANCODE_SPACE])
		{
			pause = !pause;
		}
		if(keyboard[SDL_SCANCODE_DOWN])
		{
		}
		if(keyboard[SDL_SCANCODE_UP])
		{
		}
        

        // clears the screen 
		SDL_SetRenderDrawColor( renderer, 0, 0, 0, 255 );
		SDL_RenderClear(renderer); 
		

        //Handles translating with mouse
        if(mouseonce[1]){
            coormousetranslation=0;

            //translates coordinate system
            coormousetranslation = 1;
            coordmouseclickdownx = mouseX;
            coordmouseclickdowny = mouseY;
            coordmouseoriginalbx = bx;
            coordmouseoriginalby = by;
        }
        if(mouse[1]){
            if(coormousetranslation==1){
                //translates coordinate system
                    double Sx = (coordmouseclickdownx)/kx-coordmouseoriginalbx;
                    double Sy = (coordmouseclickdowny)/ky-coordmouseoriginalby;
                    bx = (mouseX)/kx - Sx;
                    by = (mouseY)/ky - Sy;
            }
        }


        //Draws Axis and related
        SDL_SetRenderDrawColor(renderer,255,255,255,255);
        //Draws axis
        CSDL_RenderDrawRoundLine(renderer, 0, ( 0+by)*ky, SCREEN_WIDTH, (0+by)*ky, SCREEN_K*0.005/3);
        CSDL_RenderDrawRoundLine(renderer, ( 0+bx)*kx, 0, (0+bx)*kx, SCREEN_HEIGHT, SCREEN_K*0.005/3);

        //Adjusts how many points will be shown on screen
        double screeninterval = SCREEN_WIDTH/kx;
        double incrementsize = 1;
        while(screeninterval/incrementsize>7) incrementsize*=2;
        while(screeninterval/incrementsize<7) incrementsize/=2;
        double jd;
        //Draws Circles representing y-axis units
        for(jd=(int)((-by)/incrementsize+1)*incrementsize;jd>(int)((SCREEN_HEIGHT/ky-by)/incrementsize-1)*incrementsize;jd-=incrementsize){
            CSDL_RenderDrawCircle(renderer,(0+bx)*kx,(jd+by)*ky,SCREEN_K*0.01/3);
            char buffer[100];
            char precision[100];
            sprintf(precision, "%%.%df", (int)fmax(0,fmin(3,(log(incrementsize)/(-log(2))))));
            sprintf(buffer, precision,jd);
            if(jd!=0)CSDL_RenderDrawText(renderer,buffer, (-30/kx+bx)*kx, (jd+by)*ky, SCREEN_WIDTH*0.02,SCREEN_WIDTH*0.003);
        }

        screeninterval = SCREEN_WIDTH/kx;
        incrementsize = 1;
        while(screeninterval/incrementsize>10) incrementsize*=2;
        while(screeninterval/incrementsize<10) incrementsize/=2;
        //Draws Circles representing x-axis units
        for(jd=(int)((-bx)/incrementsize-1)*incrementsize;jd<(int)((SCREEN_WIDTH/kx-bx)/incrementsize+1)*incrementsize;jd+=incrementsize){
            CSDL_RenderDrawCircle(renderer,(jd+bx)*kx,(0+by)*ky,SCREEN_K*0.01/3);
            char buffer[100];
            char precision[100];
            sprintf(precision, "%%.%df", (int)fmax(0,fmin(3,(log(incrementsize)/(-log(2))))));
            sprintf(buffer, precision,jd);
            if(jd!=0)CSDL_RenderDrawText(renderer,buffer, (jd+bx)*kx, (30/ky+by)*ky, SCREEN_WIDTH*0.02,SCREEN_WIDTH*0.003);
        }

        //Runs annealing
        for(i=0;i<hc->n;i++){
            if(!pause || keyboardonce[SDL_SCANCODE_N])hc = simulatedAnnealing(hc, temp, simulatedAnnealingSuccessorCallback, simulatedAnnealingScoreCallback, simulatedAnnealingDestroyCallback);
            if(HCScore(hc) > HCScore(hcmax)){
                HCDestroy(hcmax);
                hcmax = HCCopy(hc);
            }
        }
        //if((!pause && counter%30) || keyboardonce[SDL_SCANCODE_N])hc = HCGetSuccessor(hc);
        //printf("Temp = %f\n",temp);
        //printf("Score = %f\n",HCScore(hc));

        char buffer[100];
        sprintf(buffer, "Score = %.2f, temp = %.4f, max score = %.2f",HCScore(hc),temp,HCScore(hcmax));
        CSDL_RenderDrawText(renderer,buffer, SCREEN_WIDTH/10*2, SCREEN_HEIGHT/10, SCREEN_K*0.01,SCREEN_K*0.0003);

        if(!pause)temp *= coolingRate;

        //Draws points
        for(i=0;i<npoints;i++){
            SDL_SetRenderDrawColor(renderer,0,0,0,255);
            CSDL_RenderDrawCircle(renderer, (points[i][0]+bx)*kx,(points[i][1]+by)*ky,SCREEN_K*0.01/3);
            SDL_SetRenderDrawColor(renderer,255,255,255,255);
            CSDL_RenderDrawCircle(renderer, (points[i][0]+bx)*kx,(points[i][1]+by)*ky,SCREEN_K*0.009/3);

        }

        //Draws current egdes
        if(!uncrossing){
            for(i=0;i<hc->n;i++){
                int u = hc->order[i];
                int v = hc->order[(i+1)%hc->n];

                SDL_SetRenderDrawColor(renderer,255,0,0,255);
                CSDL_RenderDrawRoundLine(renderer, (points[u][0]+bx)*kx,(points[u][1]+by)*ky, (points[v][0]+bx)*kx,(points[v][1]+by)*ky, SCREEN_K*0.01/3);
            }
        }

        //Draws max egdes
        if(!uncrossing){
            for(i=0;i<hcmax->n;i++){
                int u = hcmax->order[i];
                int v = hcmax->order[(i+1)%hcmax->n];

                SDL_SetRenderDrawColor(renderer,0,0,255,255);
                CSDL_RenderDrawRoundLine(renderer, (points[u][0]+bx)*kx,(points[u][1]+by)*ky, (points[v][0]+bx)*kx,(points[v][1]+by)*ky, SCREEN_K*0.01/3);
            }
        }

        //Performs the uncrossing operation
        if(!uncrossing && keyboardonce[SDL_SCANCODE_U]){
            uncrossing = 1;
            hcuncross = HCCopy(hcmax);
        }


        if(uncrossing){

            if(canbeuncrossed){
                canbeuncrossed=1;
                int ** mrep = (int**)malloc(hcuncross->n*sizeof(int*));
                for(i=0;i<hcuncross->n;i++){
                    mrep[i] = (int*)calloc(hcuncross->n,sizeof(int));
                }
                for(i=0;i<hcuncross->n;i++){
                    mrep[hcuncross->order[i]][hcuncross->order[(i+1)%hcuncross->n]]=1;
                    mrep[hcuncross->order[(i+1)%hcuncross->n]][hcuncross->order[i]]=1;
                }
                //printf("BEFORE: \n");
                //printf("List ");
                //for(i=0;i<hcuncross->n;i++){
                //    printf("%d ",hcuncross->order[i]);
                //}
                //printf("\nMatrix\n");
                //for(i=0;i<hcuncross->n;i++){
                //    for(j=0;j<hcuncross->n;j++){
                //        printf("%d ",mrep[i][j]);
                //    }
                //    printf("\n");
                //}
                int flag = 1;
                for(i=0;flag && i<hcuncross->n;i++){
                    int u = hcuncross->order[i];
                    int v = hcuncross->order[(i+1)%hcuncross->n];
                    for(j=0;flag && j<hcuncross->n;j++){
                        int w = hcuncross->order[j];
                        int z = hcuncross->order[(j+1)%hcuncross->n];
                        if(u!=w && u!=z && v!=w && v!=z){
                            if(weights[u][v]+weights[w][z] > weights[u][w]+weights[v][z]){
                                canbeuncrossed=1;
                                //printf("Switching (%d,%d) (%d,%d) to (%d,%d) (%d,%d)\n",u,v,w,z,u,w,v,z);
                                mrep[u][v] = 0;
                                mrep[v][u] = 0;
                                mrep[w][z] = 0;
                                mrep[z][w] = 0;

                                mrep[u][w]=1;
                                mrep[w][u]=1;
                                mrep[v][z]=1;
                                mrep[z][v]=1;

                                flag = 0;
                                break;
                            }
                        }
                    }
                }
                //recontructs list representation
                int current = 0;
                int prev = -1;
                for(i=0;i<hcuncross->n;i++){
                    hcuncross->order[i] = current;
                    for(j=0;j<hcuncross->n;j++){
                        if(j!= prev && (mrep[current][j] || mrep[j][current])){
                            prev = current;
                            current = j;
                            break;
                        }
                    }
                }

                for(i=0;i<hcuncross->n;i++){
                    free(mrep[i]);
                }
                free(mrep);

            }
            char buffer[100];
            sprintf(buffer, "uncrossed score = %.2f",HCScore(hcuncross));
            SDL_SetRenderDrawColor(renderer,255,255,255,255);
            CSDL_RenderDrawText(renderer,buffer, SCREEN_WIDTH/10, SCREEN_HEIGHT/10*1.5, SCREEN_K*0.01,SCREEN_K*0.0003);
            //Draws max egdes
            for(i=0;i<hcuncross->n;i++){
                int u = hcuncross->order[i];
                int v = hcuncross->order[(i+1)%hcuncross->n];

                SDL_SetRenderDrawColor(renderer,0,255,0,255);
                CSDL_RenderDrawRoundLine(renderer, (points[u][0]+bx)*kx,(points[u][1]+by)*ky, (points[v][0]+bx)*kx,(points[v][1]+by)*ky, SCREEN_K*0.01/3);
            }
        }

        //Example of how to draw a curve
        //Drawing a curve
        //SDL_SetRenderDrawColor(renderer,0,0,255,255);
        //double dx = SCREEN_WIDTH/kx*0.005;
        //double x = -bx-1-dx;
        ////double a = -fmax((double)(counter%200)-100.0,-(double)(counter%200)+100.0)/100.0*10;
        //double a =-4;
        //double fpx = pow(x,-a);
        //for(x=-bx-1;x<SCREEN_WIDTH/kx-bx+1;){
        //    double fx = pow(x,a);
        //    //CSDL_RenderDrawCircle(renderer,(j+bx)*kx,(j*j+by)*ky,3);
        //    if(!(fx!=fx || fpx != fpx)){
        //        if(((fx+by)*ky>0&&(fx+by)*ky<SCREEN_HEIGHT) || ((fpx+by)*ky>0&&(fpx+by)*ky<SCREEN_HEIGHT)){
        //                CSDL_RenderDrawRoundLine(renderer, ( x+bx)*kx, (fx+by)*ky, (x-dx+bx)*kx, (fpx+by)*ky, 5);
        //        }
        //    }
        //    fpx=fx;
        //    x+=dx;
        //}

        //Displays coordinates of mouse pointer
        if(mouse[2]){
            SDL_SetRenderDrawColor( renderer, 255, 255, 255, 255 );
            char buffer[100];
            sprintf(buffer, "(%.4lf,%.4lf)",(mouseX/kx-bx),(mouseY/ky-by));
            CSDL_RenderDrawText(renderer,buffer,mouseX,mouseY-30, SCREEN_WIDTH*0.02,SCREEN_WIDTH*0.003);
            CSDL_RenderDrawCircle(renderer,mouseX,mouseY,10);
        }

        //Text Display Example
//		char buffer[100]="Test";
//		SDL_SetRenderDrawColor(renderer,255,255,255,255);
//		CSDL_RenderDrawText(renderer,buffer, SCREEN_WIDTH/2, SCREEN_HEIGHT/15, SCREEN_WIDTH/30, 3);
        
		SDL_RenderPresent(renderer); 

        //Handles fps
//		calculates to 60 fps 
		counter++;
//		SDL_Delay(1000 / FPS);
	} 

    //SDL cleanup
	// destroy renderer 
	SDL_DestroyRenderer(renderer); 

	// destroy window 
	SDL_DestroyWindow(window);	

	return 0; 
} 
