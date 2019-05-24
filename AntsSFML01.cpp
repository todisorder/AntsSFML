// g++ -I./SFML/include -std=c++11 -framework sfml-window -framework sfml-graphics -framework sfml-system  AntsSFML01.cpp -o ants
// g++ -std=c++11 AntsSFML01.cpp -lsfml-window -lsfml-graphics  -lsfml-system -o ants
////////////////////////////////////////////////////////////
// Headers
////////////////////////////////////////////////////////////
#include <SFML/Graphics.hpp>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <random>

#include "resources/parameters.cpp"
#include "resources/auxiliary_functions.cpp"

#define Here(a) std::cout << "output: " #a << '\n'
#define Here2(a,b) std::cout << (a) << (b) << '\n'

void print_this(float whatever)
{
    std::cout << whatever << '\n';
}


////////////////////////////////////////////////////////////
// Class for droplets
////////////////////////////////////////////////////////////

class Droplet : public sf::Drawable
{
public:
    static int number_of_droplets;

    Droplet()
    {
    }
    Droplet(float time, float x, float y)
    {
        number_of_droplets++;
        origin_time = time;
        pos.x = x;
        pos.y = y;
        float ballRadius = 3.f;
        body.setRadius(ballRadius);
        body.setOutlineThickness(0);
        body.setOutlineColor(sf::Color::Black);
        body.setFillColor(sf::Color::White);
        body.setOrigin(ballRadius / 2, ballRadius / 2);
        body.setPosition(x,y);
        //        std::cout << number_of_droplets << "\n";
        
    }
    sf::CircleShape body;
    float origin_time;
    sf::Vector2f pos;
    float elapsed_time(float time)
    {
        return time - origin_time;
    }
    void draw(sf::RenderTarget& target, sf::RenderStates states) const
    {
        
    }
    int static CountDroplets()
    {
        return number_of_droplets;
    }
    void initialize(float time, float x, float y)   // not used...
    {
        origin_time = time;
        pos.x = x;
        pos.y = y;
        float ballRadius = 2.f;
        body.setRadius(ballRadius);
        body.setOutlineThickness(0);
        body.setOutlineColor(sf::Color::Black);
        body.setFillColor(sf::Color::White);
        body.setOrigin(ballRadius / 2, ballRadius / 2);
        body.setPosition(x,y);
    }
};

int Droplet::number_of_droplets = 0;




////////////////////////////////////////////////////////////
// Class for one ant
////////////////////////////////////////////////////////////
int what=1000;
class OneAnt : public sf::Drawable
{
public:
    
//    static int ants_that_were_born;
    
    OneAnt(int inputs, int nodes, int outputs) : Wheights1(inputs,std::vector<float>(nodes)), Wheights2(nodes,std::vector<float>(outputs))
    {
    }
    
    sf::CircleShape body;
    sf::CircleShape sensingradius;
    sf::Vector2f pos;
    sf::Vector2f pos_old;
    sf::Vector2f antenna_L;
    sf::Vector2f antenna_R;
    sf::CircleShape antenna_L_shape;
    sf::CircleShape antenna_R_shape;
    sf::Vector2f sting_pos;
    sf::Vector2f vel;
    sf::Vector2f vel_old;
    sf::Vector2f force;
    float sensing_angle;
    float stamina;
    float hitpoints;
    
    float felt_dropletsL;
    float felt_dropletsR;
    float felt_bodiesL_k1;
    float felt_bodiesR_k1;
    float felt_bodiesL_k2;
    float felt_bodiesR_k2;
    
    int kind;
    
    sf::Color kindColor;

    
    std::vector<float> InputState; // = {felt_dropletsL, felt_dropletsR, felt_bodiesL, felt_bodiesR, body.getRadius(), stamina, hitpoints, sensing_angle};
    std::vector<float> OutputState; // = {fpl, fpr}
    std::vector<float> Nodes;
    std::vector<std::vector<float>> Wheights1;
    std::vector<std::vector<float>> Wheights2;
    float bias;

    sf::Text statMessage;
    
    
    
//    std::vector<std::vector<float>> Wheights1(int number_of_inputs, std::vector<float> number_of_nodes);
//    std::vector<std::vector<float>> Wheights2(int number_of_nodes, std::vector<float> number_of_outputs);
    //  See https://stackoverflow.com/questions/12375591/vector-of-vectors-to-create-matrix
    

    
    bool I_am_alive;
    
    bool load_genesis(int I_want_the_kind_to_be=-1)
    {
        // Create an ant
        statMessage.setString("fdfd");
        bias = initial_bias;
        Nodes.resize(number_of_nodes);
        OutputState.resize(number_of_outputs);
        
        hitpoints = 100.;
        I_am_alive = true;
        sensing_angle = sensing_area_half_angle;
        std::uniform_int_distribution<int> unif(1,2);
        
        if (I_want_the_kind_to_be == -1) {
            kind = unif(generator);
        }
        else {
            kind = I_want_the_kind_to_be;
        }
        if (kind == 1) {
            stamina = 50.;
        } else {
            stamina = 1.;
        }


        
        float ballRadius = zoom_multiplier*10.f;
        body.setRadius(ballRadius - 0);
        body.setOutlineThickness(0);
        body.setOutlineColor(sf::Color::Black);
        if (kind==1) {
            kindColor = sf::Color(200,55,55);
        } else {
            kindColor = sf::Color(55,55,250);
        }
        body.setFillColor(kindColor);
        
        body.setOrigin(ballRadius / 2.  , ballRadius / 2.);
        
        float testballRadius = zoom_multiplier*30.f;
        sensingradius.setRadius(testballRadius - 0);
        sensingradius.setOutlineThickness(0);
        sensingradius.setOutlineColor(sf::Color::Black);
        sensingradius.setFillColor(sf::Color(200,255,55,50));
        sensingradius.setOrigin(testballRadius / 2. , testballRadius / 2.);
        
        float antennaRadius = zoom_multiplier*3.f;
        antenna_L_shape.setRadius(antennaRadius - 0);
        antenna_L_shape.setOutlineThickness(0);
        antenna_L_shape.setOutlineColor(sf::Color::Black);
        antenna_L_shape.setFillColor(sf::Color::Yellow);
        antenna_L_shape.setOrigin(antennaRadius / 2, antennaRadius / 2);
        
        antenna_R_shape.setRadius(antennaRadius - 0);
        antenna_R_shape.setOutlineThickness(0);
        antenna_R_shape.setOutlineColor(sf::Color::Black);
        antenna_R_shape.setFillColor(sf::Color::Yellow);
        antenna_R_shape.setOrigin(antennaRadius / 2, antennaRadius / 2);
        
        //        float x = static_cast<float>(std::rand() % 1000);       // just to see if I can define this outside the class. I can.
        //        float y = static_cast<float>(std::rand() % 1000);
        float x = x_2*Uniform(generator) + x_1*(1.-Uniform(generator));
        float y = y_2*Uniform(generator) + y_1*(1.-Uniform(generator));
        pos.x = x;
        pos.y = y;
        pos_old.x = pos.x;
        pos_old.y = pos.y;
        body.setPosition(pos.x-ballRadius / 2.,pos.y-ballRadius / 2.);
        sensingradius.setPosition(pos.x-testballRadius / 2.,pos.y-testballRadius / 2.);
        //        x = static_cast<float>(std::rand() % 1000);
        //        y = static_cast<float>(std::rand() % 1000);
        float iniangle = UniformAngle(generator);
        vel.x = cos(iniangle);
        vel.y = sin(iniangle);
        float norm = sqrt(vel.x*vel.x + vel.y*vel.y)/100.;
        vel.x /= norm;
        vel.y /= norm;
        vel_old.x = vel.x;
        vel_old.y = vel.y;
        
        float xl = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)-sensing_angle);
        float yl = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)-sensing_angle);
        float xr = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)+sensing_angle);
        float yr = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)+sensing_angle);
        antenna_L.x = xl;
        antenna_L.y = yl;
        antenna_R.x = xr;
        antenna_R.y = yr;
        antenna_L_shape.setPosition(xl,yl);
        antenna_R_shape.setPosition(xr,yr);
        
        
        for (int i=0; i<number_of_inputs; i++){
            for (int j=0; j<number_of_nodes; j++)
            {
                Wheights1[i][j] = 0.1*Uniform(generator);
            }
        }
        for (int i=0; i<number_of_nodes; i++){
            for (int j=0; j<number_of_outputs; j++)
            {
                Wheights2[i][j] = 0.1*Uniform(generator);
            }
        }
        
        InputState = {felt_dropletsL, felt_dropletsR, felt_bodiesL_k1, felt_bodiesL_k2, felt_bodiesR_k1, felt_bodiesR_k2, body.getRadius(), stamina, hitpoints, sensing_angle};
        
        return true;
    }
    bool load()
    {
        // Create an ant
        
        Nodes.resize(number_of_nodes);
        OutputState.resize(number_of_outputs);
//
        if (kind == 1) {
            stamina = 50.;
        } else {
            stamina = 1.;
        }
//        stamina = 1.;
        hitpoints = 100.;
        I_am_alive = true;
        
        float ballRadius = body.getRadius();
        body.setOutlineThickness(0);
        body.setOutlineColor(sf::Color::Black);
        if (kind==1) {
            kindColor = sf::Color(200,55,55);
        } else {
            kindColor = sf::Color(55,55,250);
        }
        body.setFillColor(kindColor);

        body.setOrigin(ballRadius / 2.  , ballRadius / 2.);
        
        float testballRadius = sensingradius.getRadius();
        sensingradius.setOutlineThickness(0);
        sensingradius.setOutlineColor(sf::Color::Black);
        sensingradius.setFillColor(sf::Color(200,255,55,30));
        sensingradius.setOrigin(testballRadius / 2. , testballRadius / 2.);
        
        float antennaRadius = zoom_multiplier*3.f;
        antenna_L_shape.setRadius(antennaRadius - 0);
        antenna_L_shape.setOutlineThickness(0);
        antenna_L_shape.setOutlineColor(sf::Color::Black);
        antenna_L_shape.setFillColor(sf::Color::Yellow);
        antenna_L_shape.setOrigin(antennaRadius / 2, antennaRadius / 2);
        
        antenna_R_shape.setRadius(antennaRadius - 0);
        antenna_R_shape.setOutlineThickness(0);
        antenna_R_shape.setOutlineColor(sf::Color::Black);
        antenna_R_shape.setFillColor(sf::Color::Yellow);
        antenna_R_shape.setOrigin(antennaRadius / 2, antennaRadius / 2);
        
//        float x = x_2*Uniform(generator) + x_1*(1.-Uniform(generator));
//        float y = y_2*Uniform(generator) + y_1*(1.-Uniform(generator));
//        pos.x = x;
//        pos.y = y;
//        pos_old.x = pos.x;
//        pos_old.y = pos.y;
        body.setPosition(pos.x-ballRadius / 2.,pos.y-ballRadius / 2.);
        sensingradius.setPosition(pos.x-testballRadius / 2.,pos.y-testballRadius / 2.);
        //        x = static_cast<float>(std::rand() % 1000);
        //        y = static_cast<float>(std::rand() % 1000);
        float iniangle = UniformAngle(generator);
        vel.x = cos(iniangle);
        vel.y = sin(iniangle);
        float norm = sqrt(vel.x*vel.x + vel.y*vel.y)/100.;
        vel.x /= norm;
        vel.y /= norm;
        vel_old.x = vel.x;
        vel_old.y = vel.y;
        
        float srad = sensingradius.getRadius();
//        std::cout <<"my bias is: " << bias << "\n";
        float xl = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)-sensing_angle);
        float yl = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)-sensing_angle);
        float xr = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)+sensing_angle);
        float yr = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)+sensing_angle);
        antenna_L.x = xl;
        antenna_L.y = yl;
        antenna_R.x = xr;
        antenna_R.y = yr;
        antenna_L_shape.setPosition(xl,yl);
        antenna_R_shape.setPosition(xr,yr);
        
        InputState = {felt_dropletsL, felt_dropletsR, felt_bodiesL_k1, felt_bodiesL_k2, felt_bodiesR_k1, felt_bodiesR_k2, body.getRadius(), stamina, hitpoints, sensing_angle};
        
        return true;
    }

    void update_ant(float time)
    {
        
        pos.x = pos_old.x + delta_t * vel_old.x;
        pos.y = pos_old.y + delta_t * vel_old.y;
        vel.x = vel_old.x + delta_t * (1./tau) * (-vel_old.x + force.x);
        vel.y = vel_old.y + delta_t * (1./tau) * (-vel_old.y + force.y);
        Boundarify();
        
        pos_old.x = pos.x;
        pos_old.y = pos.y;
        vel_old.x = vel.x;
        vel_old.y = vel.y;
        float srad = sensingradius.getRadius();
        antenna_L.x = pos.x + srad * cos(Angle(vel.x,vel.y)-sensing_angle);
        antenna_L.y = pos.y + srad * sin(Angle(vel.x,vel.y)-sensing_angle);
        antenna_R.x = pos.x + srad * cos(Angle(vel.x,vel.y)+sensing_angle);
        antenna_R.y = pos.y + srad * sin(Angle(vel.x,vel.y)+sensing_angle);
        
        float rad = body.getRadius();
        
        body.setOrigin(rad/2.,rad/2.);
        body.setPosition(pos.x-rad/2., pos.y-rad/2.);
        sf::Color tmpcol = body.getFillColor();
        tmpcol.a = 255*oneminusexp(stamina+0.4);
//        sf::Color tmpcol(0,0,0,255*oneminusexp(stamina+0.4));
        body.setFillColor(tmpcol);
//        body.setFillColor(sf::Color(200,55,55,(oneminusexp(stamina+0.4))*255));
        float radt = sensingradius.getRadius()*1.00;
        sensingradius.setRadius(radt);
        sensingradius.setOrigin(radt/2.,radt/2.);
        sensingradius.setPosition(pos.x-radt/2., pos.y-radt/2.);
        antenna_L_shape.setPosition(antenna_L.x,antenna_L.y);
        antenna_R_shape.setPosition(antenna_R.x,antenna_R.y);
        
        float distance_covered = sqrt(delta_t * vel_old.x*delta_t * vel_old.x + delta_t * vel_old.y *delta_t * vel_old.y);
        stamina = (1.-0.004*(1.+distance_covered)*(1.+0.005*body.getRadius()*body.getRadius()))*stamina;
        
//        InputState = {felt_dropletsL, felt_dropletsR, felt_bodiesL, felt_bodiesR, body.getRadius(), stamina, hitpoints, sensing_angle};

        
    }
    
    void draw(sf::RenderTarget& target, sf::RenderStates states) const
    {
        
    }
    
    float VeryComplicatedFunction();
    
    void ComputeForce(float time);
    float Feeling(float time, float x, float y, std::vector<Droplet>& pheromone);
    void Boundarify();
    
    void reproduce(std::vector<OneAnt*>& ants);
    void will_I_die();
    void kill(int antnumber, std::vector<OneAnt*>& ants);
    bool will_I_reproduce();
    
    void ShowPosition()
    {
        std::cout << "x = "<< pos.x << ", y = " << pos.y << "\n";
    }

//private:
    
};


////////////////////////////////////////////////////////////
// Class for Universe
////////////////////////////////////////////////////////////
class Universe
{
public:
    std::vector<OneAnt*> population;
    std::vector<Droplet> droplets = {};
    
    void initialize()
    {
        for (int i = 0; i<=number_of_ants; i++)
        {
            population.push_back(new OneAnt(number_of_inputs,number_of_nodes,number_of_outputs));
        }
        
        for (std::size_t i = 0; i < population.size(); ++i)
        {
            population[i]->load_genesis();
        }
    }
    
    void update_universe()
    {
        
    }
    
};




////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////

void OneAnt::ComputeForce(float time)
{
    float ax = pos.x;
    float ay = pos.y;

    float alx = antenna_L.x;
    float aly = antenna_L.y;
    float arx = antenna_R.x;
    float ary = antenna_R.y;
    
    float antenna_length = sqrt((alx - ax)*(alx - ax) + (aly - ay)*(aly - ay));
    
    float fpl,fpr;
    OutputState[0]=0.f;
    OutputState[1]=0.f;
    
    
    //  First layer
    for (int j=0; j<number_of_nodes; j++) {
        Nodes[j] = 0.f;
    }
    for (int i=0; i<number_of_inputs; i++){
        for (int j=0; j<number_of_nodes; j++)
        {
            Nodes[j] += InputState[i]*Wheights1[i][j];
//            std::cout << "w1["<<i<<"] = "<< Wheights1[i][j] <<"\n";
        }
    }
    //  Output
    for (int j=0; j<number_of_nodes; j++) {
        Nodes[j] = sigmoid(Nodes[j], bias);
//        std::cout << "node["<<j<<"] = "<< Nodes[j] <<"\n";

    }
    for (int i=0; i<number_of_nodes; i++){
        for (int j=0; j<number_of_outputs; j++)
        {
            OutputState[j] += Nodes[i]*Wheights2[i][j];
        }
    }
    

    
    
    fpl = sigmoid(OutputState[0],bias);
    
    fpr = sigmoid(OutputState[1],bias);

//    std::cout << "fpl = "<< fpl <<"\n";
//    fpl =1.f;
//    fpr =1.f;
//
    
    
//    float fpl = Feeling(time, antenna_L.x,antenna_L.y,pheromone);
//    float fpr = Feeling(time, antenna_R.x,antenna_R.y,pheromone);
    float denom = (fpl + fpr)*antenna_length*0.007;
    float numerx = (alx - ax)*fpl + (arx - ax)*fpr;
    float numery = (aly - ay)*fpl + (ary - ay)*fpr;
    force.x = numerx/denom;
    force.y = numery/denom;
    force.x = oneminusexp(10.*stamina)*negative_exponential(1.,.1,body.getRadius())*force.x;
    force.y = oneminusexp(10.*stamina)*negative_exponential(1.,.1,body.getRadius())*force.y;

}

float OneAnt::Feeling(float time, float x, float y, std::vector<Droplet>& pheromone)
{
    float result = 0.;
    float dx, dy;
    float elapsed_time;
    // neighbor_squares is defined in parameters.cpp:
    // std::vector<std::vector<float>> neighbor_squares = {N,S,W,E,NW,NE,SW,SE};
    
    int NN = Droplet::CountDroplets();
    for (int i=0; i<NN; ++i)
    {
        // Current droplet
        Droplet& dd = pheromone[i];
        // Distance from droplet to current antenna
        dx = x-dd.pos.x;
        dy = y-dd.pos.y;
        elapsed_time = time - dd.origin_time;
        if (BoundaryMethod=="periodic")
        {
            float dist = (dx*dx + dy*dy);
            if (dist < body.getRadius()*body.getRadius())
            {
                std::cout << "sdf"<<"\n";
                result += 0*Amount;
                // para comerem phero, aumenta stamina, etc... - isto vai para o main acho eu
            }
        }
        if (BoundaryMethod=="reflective")
        {
            float dist = (dx*dx + dy*dy);
            if (dist*0 <= ignore_droplet_if_this_far2)    // Careful... this does not measure across periodic Bdry!
            {
                result += Heat(dx,dy,elapsed_time,Amount);
            }
        }
        
    }
//    print_this(result);
    return std::max(threshold,result);
}

void OneAnt::reproduce(std::vector<OneAnt*>& ants)
{
    std::cout <<"reproduced havig stamina "<<stamina<<"\n";
    ants.push_back(new OneAnt(number_of_inputs,number_of_nodes,number_of_outputs));
    OneAnt* newborn = ants[ants.size()-1];
    
    newborn->kind = kind;
    newborn->pos.x = pos.x;
    newborn->pos.y = pos.y;
    newborn->pos_old.x = pos_old.x;
    newborn->pos_old.y = pos_old.y;
    newborn->body.setRadius(Perturb(body.getRadius()));
    newborn->sensingradius.setRadius(Perturb(sensingradius.getRadius()));
    newborn->sensing_angle = Perturb(sensing_angle);
    newborn->bias = PerturbWithoutSign(bias);
//    std::cout << "sa = "<< sensing_angle <<"\n";
    for (int i=0; i<number_of_inputs; i++) {
        for (int j=0; j<number_of_nodes; j++) {
            float sdfds = PerturbWithoutSign(Wheights1[i][j]);
            newborn->Wheights1[i][j] = sdfds;
//            std::cout << "just inherited w "<< Wheights1[i][j] <<" which became "<< newborn->Wheights1[i][j] <<"\n";
        }
    }
    for (int i=0; i<number_of_nodes; i++) {
        for (int j=0; j<number_of_outputs; j++) {
            float sdfds = PerturbWithoutSign(Wheights2[i][j]);
            newborn->Wheights2[i][j] = sdfds;
        }
    }
    newborn->load();
    stamina *= 0.5;
}

void OneAnt::will_I_die()
{
//    if (hitpoints <= 1. || !rolldice(std::floor(stamina)*100,100))
//    {
////        I_am_alive = false;
//    }
    
    I_am_alive = (stamina +0.05*Uniform(generator) > .005);
    
    
    // AJUSTAR ISTO:
//    I_am_alive = rolldice(std::floor(oneminusexp(stamina)*100)+1,5);
//    std::cout << std::floor(oneminusexp(stamina)*100)+1<< "\n";
    //juntar stamina
}




bool OneAnt::will_I_reproduce()
{
//    bool result = rolldice(1,700);
    bool result = rolldice(std::floor(oneminusexp(3.*stamina)*100),9000);
//    if (rolldice(std::floor(100.f-oneminusexp(stamina)*100),100)) {
//        result = false;
//    }
    return result;
}

void OneAnt::kill(int antnumber, std::vector<OneAnt*>& ants)
{
    ants.erase(ants.begin()+antnumber);
}

void OneAnt::Boundarify()
{
    if (BoundaryMethod == "reflective")
    {
        float newposx = std::max(x_1,std::min(x_2,pos.x));
        if (newposx >= x_2 || newposx <= x_1)
        {
            vel.x = -vel.x;
        }
        float newposy = std::max(y_1,std::min(y_2,pos.y));
        if (newposy >= y_2 || newposy <= y_1)
        {
            vel.y = -vel.y;
        }
    }
    if (BoundaryMethod == "periodic")
    {
        Periodify(pos.x,pos.y);
        Periodify(antenna_L.x,antenna_L.y);
        Periodify(antenna_R.x,antenna_R.y);
    }
    
}

////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////

int main()
{
    // Define some constants
    const float pi = 3.14159f;
//    float Lx = x_2-x_1;
//    float Ly = y_2-y_1;
    const int windowWidth = Lx;
    const int windowHeight = Ly;
//    float ballRadius = 10.f;
//    int number_of_ants = 35;  // in another file.
    number_of_ants -= 1;
    
    
    Universe world;
    world.initialize();
    std::vector<OneAnt*> all_the_ants = world.population;
    
    // Create the Ants
//    std::vector<OneAnt*> all_the_ants;
//    for (int i = 0; i<=number_of_ants; i++) {
//        all_the_ants.push_back(new OneAnt());
//    }
//    std::size_t current = 0;    // ?????
//
//    for (std::size_t i = 0; i < all_the_ants.size(); ++i)
//    {
//        all_the_ants[i]->load();
//    }
    // Create the pheromone vector
    std::vector<Droplet> all_droplets = world.droplets;

    
    // Create the window of the application
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight, 32), "Evol",
                            sf::Style::Titlebar | sf::Style::Close);
    window.setVerticalSyncEnabled(true);
    
    // Load the text font
    sf::Font font;
    if (!font.loadFromFile("resources/sansation.ttf"))
    return EXIT_FAILURE;
    
    // Initialize the pause message
    sf::Text pauseMessage;
    pauseMessage.setFont(font);
    pauseMessage.setCharacterSize(40);
    pauseMessage.setPosition(170.f, 150.f);
    pauseMessage.setFillColor(sf::Color::White);
    pauseMessage.setString("Press space to start/pause\n d - deposit pheromone\n r - remove pheromone\n s - save screenshot");
    
    sf::Text stats;
    stats.setFont(font);
    stats.setCharacterSize(20);
    pauseMessage.setPosition(170.f, 250.f);
    stats.setFillColor(sf::Color::White);
    


    sf::Clock clock;
    sf::Clock pause_clock;
    
    sf::Time re_start_time;
    sf::Time I_have_been_paused_this_time;
    sf::Time exec_time;
    sf::Time last_pause_time;
    
    
    bool I_will_drop = false;
    float last_drop_time = 0.;
    float last_drop_iter = 0.;
    float drop_every_sec = 1.;
//    float drop_every_iter = 10.;
    int iteration = 0;
    bool isPaused = true;
    int screenshot_count = 1;
    
    // Main loop
    while (window.isOpen())
    {
        float mouse_x = sf::Mouse::getPosition(window).x;
        float mouse_y = sf::Mouse::getPosition(window).y;

//        pauseMessage.setPosition(mx,my);
//        window.draw(pauseMessage);
        for (std::size_t i = 0; i < all_the_ants.size(); ++i)
        {
            
            float dx = all_the_ants[i]->pos.x - mouse_x;
            float dy = all_the_ants[i]->pos.y - mouse_y;
            float dist2 = (dx*dx + dy*dy);
            if (sqrt(dist2) <= all_the_ants[i]->sensingradius.getRadius())
            {
                stats.setPosition(mouse_x+20,mouse_y-10);
                stats.setString("Nr "+std::to_string(i)+"\nFDL "+std::to_string(all_the_ants[i]->felt_dropletsL)+"\nFDR "+std::to_string(all_the_ants[i]->felt_dropletsR)+"\nFBLk1 "+std::to_string(all_the_ants[i]->felt_bodiesL_k1)+"\nFBRk1 "+std::to_string(all_the_ants[i]->felt_bodiesR_k1)+"\nFBLk2 "+std::to_string(all_the_ants[i]->felt_bodiesL_k2)+"\nFBRk2 "+std::to_string(all_the_ants[i]->felt_bodiesR_k2)+"\nSTA "+std::to_string(all_the_ants[i]->stamina)+"\nOutL "+std::to_string(all_the_ants[i]->OutputState[0])+"\nOutR "+std::to_string(all_the_ants[i]->OutputState[1])+"\nBIA "+std::to_string(all_the_ants[i]->bias)+"\nHP "+std::to_string(all_the_ants[i]->hitpoints));
            }
        }
        
        
        // Handle events
        sf::Event event;
        while (window.pollEvent(event))
        {
            // Window closed or escape key pressed: exit
            if ((event.type == sf::Event::Closed) ||
                ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Escape)))
            {
                window.close();
                break;
            }
            // Space key pressed: play
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::Space))
            {
                if (isPaused)
                {
                    // (re)start the simulation
                    isPaused = false;
                    re_start_time = clock.getElapsedTime();
                    I_have_been_paused_this_time += pause_clock.getElapsedTime();
                    
//                    clock.restart();
                }
                else
                {
                    isPaused = true;
                    last_pause_time = pause_clock.getElapsedTime();
                    pause_clock.restart();
                    window.draw(pauseMessage);
                }
            }
            // Insert pheromone with mouse and keys with D
            if ((event.type == sf::Event::MouseButtonPressed) && (event.mouseButton.button == sf::Mouse::Left))
            {
//                float current_time = clock.getElapsedTime().asSeconds() - I_have_been_paused_this_time.asSeconds();
                float current_time = exec_time.asSeconds();
                all_droplets.push_back(Droplet(current_time,event.mouseButton.x,event.mouseButton.y));
            }
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::D))
            {
//                float current_time = clock.getElapsedTime().asSeconds() - I_have_been_paused_this_time.asSeconds();
                float current_time = exec_time.asSeconds();
            all_droplets.push_back(Droplet(current_time,sf::Mouse::getPosition(window).x,sf::Mouse::getPosition(window).y));
            }
            // Create new ants with f
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::F))
            {
                all_the_ants.push_back(new OneAnt(number_of_inputs,number_of_nodes,number_of_outputs));
                OneAnt* newborn = all_the_ants[all_the_ants.size()-1];
                newborn->load_genesis(1);
                newborn->pos.x = sf::Mouse::getPosition(window).x;
                newborn->pos.y = sf::Mouse::getPosition(window).y;
                newborn->pos_old.x = newborn->pos.x;
                newborn->pos_old.y = newborn->pos.y;
            }
            // Create new prey with g
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::G))
            {
                all_the_ants.push_back(new OneAnt(number_of_inputs,number_of_nodes,number_of_outputs));
                OneAnt* newborn = all_the_ants[all_the_ants.size()-1];
                newborn->load_genesis(2);
                newborn->pos.x = sf::Mouse::getPosition(window).x;
                newborn->pos.y = sf::Mouse::getPosition(window).y;
                newborn->pos_old.x = newborn->pos.x;
                newborn->pos_old.y = newborn->pos.y;
            }
            // Erase pheromone with R
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::R))
            {
                float mouse_x = sf::Mouse::getPosition(window).x;
                float mouse_y = sf::Mouse::getPosition(window).y;
                float eraser_size = 30.;
                float eraser_size2 = eraser_size*eraser_size;
                float rightnow = clock.getElapsedTime().asSeconds();
                for (int ii=0; ii<Droplet::CountDroplets(); ii++)
                {
                    Droplet& dd = all_droplets[ii];
                    float dx = dd.pos.x - mouse_x;
                    float dy = dd.pos.y - mouse_y;
                    float dist2 = (dx*dx + dy*dy);
                    
                    float elt = rightnow-all_droplets[ii].origin_time;
                    if (dist2 <= eraser_size2)
                    {
                        all_droplets.erase(all_droplets.begin()+ii);
                        Droplet::number_of_droplets--;
                    }
                }
            }
            // Save screenshot with S
            if ((event.type == sf::Event::KeyPressed) && (event.key.code == sf::Keyboard::S))
            {
                sf::Vector2u windowSize = window.getSize();
                sf::Texture texture;
                texture.create(windowSize.x, windowSize.y);
                texture.update(window);
                sf::Image screenshot = texture.copyToImage();
                if (!screenshot.saveToFile(std::to_string(screenshot_count)+"result.png"))
                {
                    return -1;
                }
                screenshot_count ++;
            }


        }
        if (!isPaused)
        {
            exec_time = clock.getElapsedTime() - I_have_been_paused_this_time;
            float rightnow = exec_time.asSeconds();
            std::cout << "it = "<<iteration <<"\n";
            std::cout << "ants = "<<all_the_ants.size() <<"\n";
            std::cout << "drops = "<<all_droplets.size() <<"\n";
            for (int i = 0; i < all_the_ants.size(); ++i)
            {
                if (iteration - last_drop_iter >= drop_every_iter && drop_every_iter >= 0)
                {
                    I_will_drop = true;
                }
                
                all_the_ants[i]->will_I_die();
                if (all_the_ants[i]->I_am_alive)
                {
                    all_the_ants[i]->felt_dropletsL = 0.f;
                    all_the_ants[i]->felt_dropletsR = 0.f;
                    all_the_ants[i]->felt_bodiesL_k1 = 0.f;
                    all_the_ants[i]->felt_bodiesR_k1 = 0.f;
                    all_the_ants[i]->felt_bodiesL_k2 = 0.f;
                    all_the_ants[i]->felt_bodiesR_k2 = 0.f;

                    // Determine inputs for neural network and compute interaction effects
                    for (int j = 0; j < all_the_ants.size(); ++j)
                    {
                        
                        float dx,dy,distsq, dist, value;
                        if (i!=j) {
                            // Detect other ants
                            dx = all_the_ants[i]->antenna_L.x - all_the_ants[j]->pos.x;
                            dy = all_the_ants[i]->antenna_L.y - all_the_ants[j]->pos.y;
                            distsq = dx*dx + dy*dy;
                            value = negative_exponential(400.,.1,sqrt(distsq));
                            //tirar o decrescimento de stamina qd reproduz
                            
                            if (all_the_ants[j]->kind == 1) {
                                all_the_ants[i]->felt_bodiesL_k1 += value;
                            }
                            if (all_the_ants[j]->kind == 2) {
                                all_the_ants[i]->felt_bodiesL_k2 += value;
                            }

//                            all_the_ants[i]->felt_bodiesL += negative_exponential(400.,.1,sqrt(distsq));
//                            std::cout <<"am adding ant "<<j<<" to feeling of ant "<<i<<"\n";
                            
                            dx = all_the_ants[i]->antenna_R.x - all_the_ants[j]->pos.x;
                            dy = all_the_ants[i]->antenna_R.y - all_the_ants[j]->pos.y;
                            distsq = dx*dx + dy*dy;
                            value = negative_exponential(400.,.1,sqrt(distsq));
                            
                            if (all_the_ants[j]->kind == 1) {
                                all_the_ants[i]->felt_bodiesR_k1 += value;
                            }
                            if (all_the_ants[j]->kind == 2) {
                                all_the_ants[i]->felt_bodiesR_k2 += value;
                            }

//                            all_the_ants[i]->felt_bodiesR += negative_exponential(400.,.1,sqrt(distsq));
                        }
                        dx = all_the_ants[i]->pos.x - all_the_ants[j]->pos.x;
                        dy = all_the_ants[i]->pos.y - all_the_ants[j]->pos.y;
                        dist = sqrt(dx*dx + dy*dy);
                        if (dist<0.5*(all_the_ants[i]->body.getRadius()+all_the_ants[j]->body.getRadius()) && i!=j )
                        {
                            if (all_the_ants[i]->kind == 1 && all_the_ants[j]->kind == 2) {
                                std::cout <<"ant "<<i<<" is eating ant "<<j<<"\n";
                                all_the_ants[i]->stamina *= 3.5;
                                all_the_ants[j]->stamina *= 0.5;
                            }
                            if (all_the_ants[i]->kind ==  all_the_ants[j]->kind) {
                                std::cout <<"ant "<<i<<" and "<<j<<" are hurting\n";
                                if (all_the_ants[i]->kind == 1) {
                                    all_the_ants[i]->stamina *= 0.3;
                                } else {
                                    all_the_ants[i]->stamina *= 0.8;
                                }
                                
//                                all_the_ants[j]->stamina *= 0.3;
                            }
                            
                        }
                        // outras coisas...
                    }
                    for (std::size_t ii=0; ii<all_droplets.size(); ii++)
                    {
                        float dx,dy,distsq;
                        dx = all_the_ants[i]->antenna_L.x - all_droplets[ii].pos.x;
                        dy = all_the_ants[i]->antenna_L.y - all_droplets[ii].pos.y;
                        distsq = dx*dx + dy*dy;
                        all_the_ants[i]->felt_dropletsL += negative_exponential(100.,.1,sqrt(distsq));
//                        all_the_ants[i]->felt_dropletsL += ((distsq < all_the_ants[i]->body.getRadius()*all_the_ants[i]->body.getRadius()*2.) ? 1.f : 0.f);

                        dx = all_the_ants[i]->antenna_R.x - all_droplets[ii].pos.x;
                        dy = all_the_ants[i]->antenna_R.y - all_droplets[ii].pos.y;
                        distsq = dx*dx + dy*dy;
                        all_the_ants[i]->felt_dropletsR += negative_exponential(100.,.1,sqrt(distsq));
//                        all_the_ants[i]->felt_dropletsR += ((distsq < all_the_ants[i]->body.getRadius()*all_the_ants[i]->body.getRadius()*2.) ? 10.f : 0.f);

//                        std::cout << "AA["<<i<<"] = "<< all_the_ants[i]->felt_dropletsR<<"\n";
                    }
                    all_the_ants[i]->InputState  = {all_the_ants[i]->felt_dropletsL, all_the_ants[i]->felt_dropletsR, all_the_ants[i]->felt_bodiesL_k1, all_the_ants[i]->felt_bodiesL_k2, all_the_ants[i]->felt_bodiesR_k1, all_the_ants[i]->felt_bodiesR_k2, all_the_ants[i]->body.getRadius(), all_the_ants[i]->stamina, all_the_ants[i]->hitpoints, all_the_ants[i]->sensing_angle};
                    
//                    std::cout << "iiiii = " << all_the_ants[i]->InputState[2];
                    
                    
                    
                
                    all_the_ants[i]->ComputeForce(rightnow);
                    all_the_ants[i]->update_ant(rightnow);
                }
                if (!all_the_ants[i]->I_am_alive)
                {
                    all_the_ants[i]->kill(i,all_the_ants);
                }
                
                if (all_the_ants[i]->will_I_reproduce())
                {
                    all_the_ants[i]->reproduce(all_the_ants);
                }
                
                for (std::size_t ii=0; ii<all_droplets.size(); ii++)
                {
                    //  Eat droplets
                    float dx = all_the_ants[i]->pos.x - all_droplets[ii].pos.x;
                    float dy = all_the_ants[i]->pos.y - all_droplets[ii].pos.y;
                    float dist = (dx*dx + dy*dy);
                    float rsq = all_the_ants[i]->body.getRadius()*all_the_ants[i]->body.getRadius();
                    if (dist < rsq && all_the_ants[i]->kind == 2)
                    {
                        all_droplets.erase(all_droplets.begin()+ii);
                        std::cout <<"Just ate drop nr "<< ii <<"\n";
                        Droplet::number_of_droplets--;
                        all_the_ants[i]->stamina *= 2.5;
                        all_the_ants[i]->hitpoints *= 1.01;
//                        float rad = all_the_ants[i]->body.getRadius();
//                        all_the_ants[i]->body.setRadius(rad*1.01);
                    }
                    //  Remove old droplets
                    float current_time = exec_time.asSeconds();
                    float elt = current_time - all_droplets[ii].origin_time;
                    if (elt > 50.) {
                        all_droplets.erase(all_droplets.begin()+ii);
                        std::cout <<"removed ph "<< ii <<"\n";
                        Droplet::number_of_droplets--;
                    }
                }
                
                if (I_will_drop)
                {
                    float current_time = exec_time.asSeconds();
                    for (int ii=0; ii<Droplet::CountDroplets(); ii++)
                    {
                        
                        float elt = current_time - all_droplets[ii].origin_time;
                        if (exp(-Evaporation * elt)*Amount <= .01)
                        {
                            
                            all_droplets.erase(all_droplets.begin()+ii);
                            std::cout <<"removed ph "<< ii <<"\n";
                            Droplet::number_of_droplets--;
                        }
                    }
                    // There is just a global pheromone:
                    all_droplets.push_back(Droplet(current_time,all_the_ants[i]->pos.x,all_the_ants[i]->pos.y));
                    last_drop_iter = iteration;
                }
            }
            // End of Ants dropping pheromone
            // Generate new droplets:
            exec_time = clock.getElapsedTime() - I_have_been_paused_this_time;
            rightnow = exec_time.asSeconds();
            bool dice = rolldice(10,14);  // boolean
            if (dice==true)
            {
                float x = x_2*Uniform(generator) + x_1*(1.-Uniform(generator));
                float y = y_2*Uniform(generator) + y_1*(1.-Uniform(generator));
                all_droplets.push_back(Droplet(rightnow,x/1.,y/1.));
                last_drop_iter = iteration;
            }
            
            
            
            
            iteration++;
            I_will_drop = false;
        }
        // Clear the window
        window.clear(sf::Color(20, 20, 20));
        
        //  Drawing
        
        // Draw the pheromone
        int N = Droplet::CountDroplets();
        //            print_this(N-all_droplets.size());
        
        for (int i=0; i<N; ++i)
        {
            // Current droplet
            Droplet& dd = all_droplets[i];
            window.draw(dd.body);
            if (!isPaused)
            {
                dd.body.setFillColor(sf::Color(255,255,255,255*exp(.1*(dd.origin_time - exec_time.asSeconds()))));
            }
        }
        
//        std::cout << "ants = " << all_the_ants.size() <<"\n";
        
        // Draw the ants
        for (std::size_t i = 0; i < all_the_ants.size(); ++i)
        {
            window.draw(all_the_ants[i]->sensingradius);
            window.draw(all_the_ants[i]->body);
            window.draw(all_the_ants[i]->antenna_L_shape);
            window.draw(all_the_ants[i]->antenna_R_shape);
        }
        if (isPaused)
        {
            window.draw(pauseMessage);
            window.draw(stats);
            
        }
        
        // Display things on screen
        window.display();
//        print_this(iteration);
    }

    return EXIT_SUCCESS;

}







