// g++ -I./SFML/include -std=c++11 -framework sfml-window -framework sfml-graphics -framework sfml-system  AntsSFML01.cpp -o ants
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
    
    OneAnt()
    {
    }
    
    sf::CircleShape body;
    sf::CircleShape testbody;
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
    
    
    bool load()
    {
        // Create an ant
        float ballRadius = zoom_multiplier*10.f;
        body.setRadius(ballRadius - 0);
        body.setOutlineThickness(0);
        body.setOutlineColor(sf::Color::Black);
        body.setFillColor(sf::Color(200,55,55));
        body.setOrigin(ballRadius / 2.  , ballRadius / 2.);

        float testballRadius = zoom_multiplier*30.f;
        testbody.setRadius(testballRadius - 0);
        testbody.setOutlineThickness(0);
        testbody.setOutlineColor(sf::Color::Black);
        testbody.setFillColor(sf::Color(200,255,55,0));
        testbody.setOrigin(testballRadius / 2. , testballRadius / 2.);

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
        testbody.setPosition(pos.x-testballRadius / 2.,pos.y-testballRadius / 2.);
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
        
        float xl = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)-sensing_area_half_angle);
        float yl = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)-sensing_area_half_angle);
        float xr = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)+sensing_area_half_angle);
        float yr = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)+sensing_area_half_angle);
        antenna_L.x = xl;
        antenna_L.y = yl;
        antenna_R.x = xr;
        antenna_R.y = yr;
        antenna_L_shape.setPosition(xl,yl);
        antenna_R_shape.setPosition(xr,yr);

        return true;
    }
    
    void update(float time, std::vector<Droplet>& pheromone)
    {
        ComputeForce(time, pheromone);
        pos.x = pos_old.x + delta_t * vel_old.x;
        pos.y = pos_old.y + delta_t * vel_old.y;
        vel.x = vel_old.x + delta_t * (1./tau) * (-vel_old.x + force.x);
        vel.y = vel_old.y + delta_t * (1./tau) * (-vel_old.y + force.y);
        Boundarify();
        pos_old.x = pos.x;
        pos_old.y = pos.y;
        vel_old.x = vel.x;
        vel_old.y = vel.y;
        antenna_L.x = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)-sensing_area_half_angle);
        antenna_L.y = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)-sensing_area_half_angle);
        antenna_R.x = pos.x + sensing_area_radius * cos(Angle(vel.x,vel.y)+sensing_area_half_angle);
        antenna_R.y = pos.y + sensing_area_radius * sin(Angle(vel.x,vel.y)+sensing_area_half_angle);
        
        float rad = body.getRadius();
        body.setPosition(pos.x-rad/2., pos.y-rad/2.);
        antenna_L_shape.setPosition(antenna_L.x,antenna_L.y);
        antenna_R_shape.setPosition(antenna_R.x,antenna_R.y);
//        body.move(vel.x,vel.y);
//        antenna_L_shape.move(vel.x,vel.y);
//        antenna_R_shape.move(vel.x,vel.y);

    }
    
    void draw(sf::RenderTarget& target, sf::RenderStates states) const
    {
        
    }
    void ComputeForce(float time, std::vector<Droplet>& pheromone);
    float FeltPheromone(float time, float x, float y, std::vector<Droplet>& pheromone);
    void Boundarify();
    void ShowPosition()
    {
        std::cout << "x = "<< pos.x << ", y = " << pos.y << "\n";
    }

//private:
    
};



////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////

void OneAnt::ComputeForce(float time, std::vector<Droplet>& pheromone)
{
//    force.x = 10.f*sin(.5*pos.x)+0.4f;
//    force.y = -15.f*cos(.5*pos.y)+0.5f;
//    force.x = 0.f;
//    force.y = 0.f;
    float ax = pos.x;
    float ay = pos.y;
    float alx = antenna_L.x;
    float aly = antenna_L.y;
    float arx = antenna_R.x;
    float ary = antenna_R.y;
    float fpl = FeltPheromone(time, antenna_L.x,antenna_L.y,pheromone);
    float fpr = FeltPheromone(time, antenna_R.x,antenna_R.y,pheromone);
    float denom = fpl + fpr;
    float numerx = (alx - ax)*fpl + (arx - ax)*fpr;
    float numery = (aly - ay)*fpl + (ary - ay)*fpr;
    force.x = numerx/denom;
    force.y = numery/denom;

}

float OneAnt::FeltPheromone(float time, float x, float y, std::vector<Droplet>& pheromone)
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
        dx = x-dd.pos.x;
        dy = y-dd.pos.y;
        elapsed_time = time - dd.origin_time;
        if (BoundaryMethod=="periodic")
        {
            float dist = (dx*dx + dy*dy);
            if (dist*0 <= ignore_droplet_if_this_far2)    // Careful... this does not measure across periodic Bdry!
            {
                result += Heat(dx,dy,elapsed_time,Amount);
                for (int dir=0; dir<neighbor_squares.size(); dir++)
                {
                    result += Heat(dx+neighbor_squares[dir][0],dy+neighbor_squares[dir][1],elapsed_time,Amount);
                }
                
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
    
    

    
    // Create the Ants
    std::vector<OneAnt*> all_the_ants;
    for (int i = 0; i<=number_of_ants; i++) {
        all_the_ants.push_back(new OneAnt());
    }
    std::size_t current = 0;    // ?????

    for (std::size_t i = 0; i < all_the_ants.size(); ++i)
    {
        all_the_ants[i]->load();
//        std::cout << all_the_ants[i]->pos.y << "\t";
    }

    // Create the pheromone vector
//    std::vector<float> droplet_centers_x;
//    std::vector<float> droplet_centers_y;
//    std::vector<float> droplet_origin_times;

    
    // Create the window of the application
    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight, 32), "SFML Ants",
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

    sf::Clock clock;
    sf::Clock pause_clock;
    
    sf::Time re_start_time;
    sf::Time I_have_been_paused_this_time;
    sf::Time exec_time;
    sf::Time last_pause_time;
    
    // I need -std=c++11 for the following:
    std::vector<Droplet> all_droplets = {};
    
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
//            std::cout << "this = " <<clock.getElapsedTime().asSeconds() - last_drop_time << "\n";
//            float rightnow = clock.getElapsedTime().asSeconds();
            exec_time = clock.getElapsedTime() - I_have_been_paused_this_time;
//            std::cout << iteration << " iterations in " << exec_time.asSeconds() << "sec." << "\r";
            float rightnow = exec_time.asSeconds();
            std::cout << all_droplets.size() << "\t" <<Droplet::CountDroplets() <<"\n";
            for (std::size_t i = 0; i < all_the_ants.size(); ++i)
            {
                if (iteration - last_drop_iter >= drop_every_iter)
                {
                    I_will_drop = true;
                }
//                all_the_ants[i]->ShowPosition();
                all_the_ants[i]->update(rightnow, all_droplets);
                // Drop pheromone
                if (I_will_drop)
                {
                    float current_time = exec_time.asSeconds();
                    for (int ii=0; ii<Droplet::CountDroplets(); ii++)
                    {
                        // Erase some droplets
//                        float elt = rightnow-all_droplets[ii].origin_time;//not used
//                        float current_time = rightnow - I_have_been_paused_this_time.asSeconds();
                     
                        
                        
                        
//                        
//                        if (exp(-Evaporation * current_time)*Amount <= .01 || current_time - all_droplets[ii].origin_time > droplet_too_old)
//                        {
//                            all_droplets.erase(all_droplets.begin()+ii);
//                            std::cout <<"removed ph"<<"\n";
//                            Droplet::number_of_droplets--;
//                        }
                        
                        
                        if (exp(-Evaporation * current_time)*Amount <= .01)
                        {
//                            std::cout << all_droplets.size() << "\t" <<Droplet::CountDroplets() <<"\n";
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            // PROBLEM HERE: all_droplets.erase(all_droplets.begin()++ii);
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            all_droplets.erase(all_droplets.begin()++ii);
                            std::cout <<"removed ph "<< ii <<"\n";
                            Droplet::number_of_droplets--;
                        }

                        
                        
                        
                        
                        
                        
                    }
                    // There is just a global pheromone:
//                    rightnow = clock.getElapsedTime().asSeconds();
//                    float current_time = rightnow - I_have_been_paused_this_time.asSeconds();
                    all_droplets.push_back(Droplet(current_time,all_the_ants[i]->pos.x,all_the_ants[i]->pos.y));
//                    last_drop_time = current_time;
                    last_drop_iter = iteration;
                }
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

        
        
        // Draw the ants
        for (std::size_t i = 0; i < all_the_ants.size(); ++i)
        {
            window.draw(all_the_ants[i]->body);
            window.draw(all_the_ants[i]->testbody);
            window.draw(all_the_ants[i]->antenna_L_shape);
            window.draw(all_the_ants[i]->antenna_R_shape);
        }
        if (isPaused)
        {
            window.draw(pauseMessage);
        }
//        // Draw the pheromone
//        int N = Droplet::CountDroplets();
////            print_this(N-all_droplets.size());
//        
//        for (int i=0; i<N; ++i)
//        {
//            // Current droplet
//            Droplet& dd = all_droplets[i];
//            window.draw(dd.body);
//            if (!isPaused)
//            {
//                dd.body.setFillColor(sf::Color(255,255,255,255*exp(.05*(dd.origin_time - clock.getElapsedTime().asSeconds()))));
//            }
//        }
        


//        else
//        {
//            // Draw the pause message
//            window.draw(pauseMessage);
//        }
        
        // Display things on screen
        window.display();
//        print_this(iteration);
    }

    return EXIT_SUCCESS;

}







