
//  Randoms
std::default_random_engine generator(seed1);
std::normal_distribution<float> Normal(0.,1.);          // Normal(0.,1.std::)
std::normal_distribution<float> SmallNormal(0.,.05);      // (0.,.05)
std::uniform_real_distribution<float> UniformAngle(0.,2.*Pi);      // Uniformly distributed angle
std::uniform_real_distribution<float> Uniform(0.,1.);      // Uniformly distributed in [0,1]
std::uniform_int_distribution<int> UniformInteger(0,10);      // Uniformly distributed integer   20
//http://www.cplusplus.com/reference/random/normal_distribution/



/////////////////////////////////////////////////
/////////////////////////////////////////////////
//      Theta de vetor  http://en.cppreference.com/w/cpp/numeric/math/atan2
//      "To compute the value, the function takes into account the sign
//      of both arguments in order to determine the quadrant."
//      Cuidado que atan2 está entre -Pi e Pi,
//      mas acho que isso não tem influencia porque
//      eu só calculo senos e cosenos do angulo,
//      que dariam a mesma coisa se fosse em (0, 2Pi).
//      CUIDADO Usage: atan2(Y,X) = arctan(Y/X) !!!!
/////////////////////////////////////////////////
/////////////////////////////////////////////////
float Angle(float X, float Y)
{
    float aux =  atan2(Y,X);
    return aux;
    
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
//      End Theta de vetor
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
//      Angle between 2 vectors (x1,y1) and (x2,y2)
//      https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
/////////////////////////////////////////////////
/////////////////////////////////////////////////
float AngleBetween(float x1, float y1, float x2, float y2)
{
    
    float dot = x1*x2 + y1*y2;
    float det = x1*y2 - y1*x2;
    float aux =  atan2(det,dot);
    return aux;
    
}


/////////////////////////////////////////////////
/////////////////////////////////////////////////
//      Heat Equation fundamental solution
/////////////////////////////////////////////////
/////////////////////////////////////////////////
float Heat(float XX, float YY, float elapsed_time, float amount)
{
    
    float aux = 0.;
    //    cout <<"????? = " << elapsed_time << endl;
    
    aux = (1. / (4.*Pi* Diffusion * elapsed_time));
    aux *= exp(-(XX*XX + YY*YY) / (4.*Diffusion*elapsed_time));
    aux *= exp(-Evaporation*elapsed_time); // Evaporation .001
    aux *= amount;
    
    return aux;
}

//  cf. http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
double Sign(double aa){
    if (aa > 0.) return 1.;
    if (aa < 0.) return -1.;
    return 0.;
}


void Periodify(float& x, float& y)
{
    if (x>=x_2 || x<=x_1)
    {
        x = x + Sign(x_2 - x)*(x_2 - x_1);
    }
    if (y>=y_2 || y<=y_1)
    {
        y = y + Sign(y_2 - y)*(y_2 - y_1);
    }
}

