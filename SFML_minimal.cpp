////////////////////////////////////////////////////////////
// Headers
////////////////////////////////////////////////////////////
#include <SFML/Graphics.hpp>
//#include <SFML/Audio.hpp>
#include <cmath>
#include <ctime>
#include <cstdlib>





int main()
{
    // Define some constants
    const float pi = 3.14159f;
    const int windowWidth = 1000;
    const int windowHeight = 1000;
    float ballRadius = 10.f;

    // Create an ant
    sf::CircleShape ball;
    ball.setRadius(ballRadius - 3);
    ball.setOutlineThickness(3);
    ball.setOutlineColor(sf::Color::Black);
    ball.setFillColor(sf::Color::White);
    ball.setOrigin(ballRadius / 2, ballRadius / 2);

    
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
    pauseMessage.setString("SFML Ants.\nPress space to start.");

    sf::Clock clock;
    bool isPlaying = false;
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
                if (!isPlaying)
                {
                    // (re)start the game
                    isPlaying = true;
                    clock.restart();
                    
                }
            }

        }
        // Clear the window
        window.clear(sf::Color(50, 100, 50));
        
        if (isPlaying)
        {
            // Draw the ball
            window.draw(ball);
        }
        else
        {
            // Draw the pause message
            window.draw(pauseMessage);
        }
        
        // Display things on screen
        window.display();
    }

    return EXIT_SUCCESS;

}