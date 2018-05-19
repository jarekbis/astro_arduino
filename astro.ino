//strefa czasowa "zimowa" czas letni jest wyliczany pozniej
#define TZ 1
#define LONGITUDE 20
#define LATITUDE 50

//for Arduino define ARDU
#define ARDU
//#undef ARDU

#include <math.h>

#ifdef ARDU

#include <avr/pgmspace.h>
#define read_float(p, o) pgm_read_float(oe+12*(p)+(o))

#else
#include <iostream>

typedef char byte;

using namespace std;
#define PROGMEM
#define PI M_PI
#define DEG_TO_RAD (PI/180)
#define read_float(p, o) (*(oe+12*(p)+(o)))
#endif // ARDU

//http://www.stjarnhimlen.se/comp/ppcomp.html

//funkcje oe2radec

float RA, Dec,
      Azimuth, Altitude, Ls, UT_body_in_south,
      Rise, Set;

#define OE_N(cbody, day) (read_float(cbody ,0)+(day)*(read_float(cbody, 1)))
#define OE_I(cbody, day) (read_float(cbody ,2)+(day)*(read_float(cbody ,3)))
#define OE_W(cbody, day) (read_float(cbody ,4)+(day)*(read_float(cbody ,5)))
#define OE_A(cbody, day) (read_float(cbody ,6)+(day)*(read_float(cbody ,7)))
#define OE_E(cbody, day) (read_float(cbody ,8)+(day)*(read_float(cbody ,9)))
#define OE_M(cbody, day) (read_float(cbody ,10)+(day)*(read_float(cbody ,11)))
//orbital elements
const PROGMEM  float oe[]  =
{
    // N(day)=N1+N2(day) longitude of ascending node -deg
    // i(day)=i1+i2(day) inclination to the ecliptic (plane of the Earth's orbit)-deg
    // w(day)=w1+w2(day) argument of perihelion -deg
    // a(day)=a1+a2(day) semi-major axis, or mean distance from Sun -AU (for moon -Earth radii)
    // e(day)=e1+e2(day) eccentricity (0 -cirle, (0-1) elipse, 1 -parabola)
    // M(day)=M1+M2(day) mean anomaly (0 at perihelion; increases uniformly with time)
    //{ N1,              N2,      i1,       i2,       w1,         w2,
    //  a1,              a2,      e1,       e2,       M1,         M2
    //},
    //Sun orbital elements
    0.0,          0.0,     0.0,      0.0, 282.9404, 4.70935E-5,
    1.000000,          0.0,0.016709,-1.151E-9, 356.0470, 0.9856002585,
    //Mercury orbital elements
    48.3313, 3.24587e-5,   7.0047, 5.00e-8,  29.1241,   1.01444e-5,
    0.387098,        0.0, 0.205635,5.59e-10, 168.6562, 4.0923344368,
    //Venus orbital elements
    76.6799, 2.46590e-5,   3.3946,   2.75e-8, 54.8910,  1.38374e-5,
    0.723330,        0.0, 0.006773, -1.302e-9, 48.0052, 1.6021302244,
    //Moon orbital elements
    125.1228, -0.0529538083, 5.1454,      0.0,318.0634, 0.1643573223,
    60.2666,           0.0,0.054900,     0.0,115.3654,13.0649929509,
    //Mars orbital elements
    49.5575,  2.11081e-5,   1.8497,  -1.78e-8,  286.5016, 2.92961e-5,
    1.523688,         0.0, 0.093405,  2.516e-9,   18.6021, 0.5240207766,
    //Juptier orbital elements
    100.4542, 2.76854e-5,   1.3030, -1.557e-7,  273.8777, 1.64505e-5,
    5.20256,        0.0, 0.048498,  4.469e-9,   19.8950, 0.0830853001,
    //Saturn orbital elements
    113.6634, 2.38980E-5,   2.4886, -1.081E-7,  339.3939, 2.97661E-5,
    9.55475,        0.0, 0.055546, -9.499E-9,  316.9670, 0.0334442282,
    //Uranus orbital elements
    74.0005,  1.3978E-5,   0.7733,    1.9E-8,   96.6612,  3.0565E-5,
    19.18171,   -1.55E-8, 0.047318,   7.45E-9,  142.5905, 0.011725806,
    //Neptune orbital elements
    131.7806,  3.0173E-5,   1.7700,  -2.55E-7,  272.8461,  -6.027E-6,
    30.05826,   3.313E-8, 0.008606,   2.15E-9,  260.2471, 0.005995147
};

double dayNumber(int year, byte m, byte D, double UT)//oblicza numer dnia od 1.1.2000r
{
    return (double)(367 * (long)year - 7 * (year + (m + 9) / 12) / 4 + 275 * m / 9 + D - 730530) + UT / 24;
}
float DST(int y, int m, int d)
{
//Czas jest zmieniany w ostatnia niedziele marca i pazdziernika odpowiednio 1:00 i 2:00 (o godz 0:00 UT)
//wyliczam tylko daty a wiec od polnocy lokalnej do godziny zmiany czasu jest "¿le" wyliczony DST
    if ((m>3)&&(m<10)) return 1;
    if (m==3)
    {
        if ( 31-((5+(int)dayNumber(y,3,31,12))%7)>d) return 0;
        else return 1;
    }
    if ((m==10) && (31-((5+(int)dayNumber(y,10,31,12))%7)>d)) return 1;

    return 0;
}
float mod360(float x)
{
    if (x>360)
        while (x >= 360)x -= 360;
    else
        while (x < 0) x += 360;
    return x;
}
void oe2radec(byte cbody, int year, byte m, byte D, double UT = 12)
{
//from orbital elements to Right Ascension (RA) and Declination (Dec)
    float d = dayNumber(year, m, D, UT);
    float M = OE_M(cbody, d);
    M = mod360(M)*DEG_TO_RAD;
    float e = OE_E(cbody, d);
    float E = M + e * sin(M) * (1.0 + e * cos(M));
    /* For comet orbits with eccentricites close to one, a difference of less than 1E-4 or 1E-5 degrees should be required.
    If this iteration formula won't converge, the eccentricity is probably too close to one. 
    Then you should instead use the formulae for near-parabolic or parabolic orbits.
    http://www.stjarnhimlen.se/comp/ppcomp.html */
    float E1;
    do {
      E1 = E;
      E = E - (E - e*sin(E) - M)/(1-e*cos(E));
    } while (abs(E1 - E) > 1e-4);
    float a = OE_A(cbody, d);
    float xv = a*(cos(E) - e);
    float yv = a*sqrt(1.0 - e * e) * sin(E);
    float v = atan2(yv, xv);
    float r = sqrt(xv * xv + yv * yv);
    float N=OE_N(cbody, d);
    N = mod360(N)*DEG_TO_RAD;
    float w=OE_W(cbody, d);
    w = mod360(w)*DEG_TO_RAD;
    float i=OE_I(cbody, d)*DEG_TO_RAD;
    float xh=r*(cos(N)*cos(v+w)-sin(N)*sin(v+w)*cos(i));
    float yh=r*(sin(N)*cos(v+w)+cos(N)*sin(v+w)*cos(i));
    float zh=r*(sin(v+w)*sin(i));
    float lonecl=atan2(yh,xh);
    float latecl=atan2(zh,sqrt(xh*xh+yh*yh));//for Sun = 0
    if (cbody == 3)
    {
        float Ls=(mod360(OE_M(0, d))-mod360(OE_W(0, d)))*DEG_TO_RAD;
        float Lm=M+N+w;
        lonecl+=(-1.274*sin(M-2*(Lm-Ls))+.658*sin(2*(Lm-Ls))-0.186*sin(mod360(OE_M(0, d))*DEG_TO_RAD))*DEG_TO_RAD;
        latecl+=(-0.173*sin(M+w-2*(Lm-Ls))+0.055*sin(w+2*(Lm-Ls)))*DEG_TO_RAD;//tutaj przy 0.055 troche poskracalem, bo dalo sie :)
        r-=(0.58*cos(M-2*(Lm-Ls))+0.46*cos(2*(Lm-Ls)));
    }
    float rs=0;
    if ((cbody !=0 ) && (cbody != 3)){
        float M = OE_M(0, d);
        M = mod360(M)*DEG_TO_RAD;
        float e = OE_E(0, d);
        float E = M + e * sin(M) * (1.0 + e * cos(M));
        rs = sqrt((cos(E) - e)*(cos(E) - e)+(1.-e*e)*sin(E)*sin(E));
    }
    Ls = (mod360(OE_M(0, d))+mod360(OE_W(0, d)))*DEG_TO_RAD;//use in radec2azimuth()
    xh=r*cos(lonecl)*cos(latecl)+rs*cos(Ls);
    yh=r*sin(lonecl)*cos(latecl)+rs*sin(Ls);
    zh=r*sin(latecl);
    float ecl= (23.4393 - 3.563E-7 * d)*DEG_TO_RAD;
    float xe=xh;
    float ye=yh*cos(ecl)-zh*sin(ecl);
    float ze=yh*sin(ecl)+zh*cos(ecl);
    RA=atan2(ye,xe);
    Dec = atan2(ze, sqrt(xe*xe+ye*ye));
}

void radec2azimuth(float UT, float longit=LONGITUDE, float lat = LATITUDE)
{
    longit *= DEG_TO_RAD; //wspolrzedne geograficzne +20E,+50N
    lat *= DEG_TO_RAD;
    UT_body_in_south = (RA - (PI + Ls) - longit) / (15 * DEG_TO_RAD);
    while (UT_body_in_south < 0) UT_body_in_south += 24;
    float HA = -UT_body_in_south+UT;
    HA*=15*DEG_TO_RAD;
    Azimuth=180+atan2(sin(HA)*cos(Dec),cos(HA)*cos(Dec)*sin(lat)-sin(Dec)*cos(lat))/DEG_TO_RAD;
    Altitude=asin(cos(HA)*cos(Dec)*cos(lat)+sin(Dec)*sin(lat))/DEG_TO_RAD;
}

void radec2riseset(float h=0, float longit=LONGITUDE, float lat=LATITUDE)
{
    h*=DEG_TO_RAD;
    longit*=DEG_TO_RAD;
    lat*=DEG_TO_RAD;
    float clha = (sin(h) - sin(lat) * sin(Dec)) / (cos(lat) * cos(Dec));
    if (abs(clha) <= 1) clha = acos(clha);
    else
    {
        clha=0;
        /*Yet another thing to consider: the Sun is always in the south near 12:00 local time, but the Moon may be in the
         * south at any time of the day (or night). This means you must pay more attention that you're really iterating
         * towards the rise or set time you want, and not some rise/set time one day earlier, or later.

        Since the Moon rises and sets on the average 50 minutes later each day, there usually will be one day each month when
        the Moon never rises, and another day when it never sets. You must have this in mind when iterating your rise/set times,
        otherwise your program may easily get caught into an infinite loop when it tries to force e.g. a rise time between
        00:00 and 24:00 local time on a day when the Moon never rises.

        At high latitudes the Moon occasionally rises, or sets, twice on a single calendar day. This may happen above the
        "lunar arctic circle", which moves between 61.5 and 71.9 deg latitude during the 18-year period of the motion of
        the lunar nodes. You may want to pay attention to this.*/
    }
    Rise = UT_body_in_south - clha * 12 / PI;
    Set = UT_body_in_south + clha * 12 / PI;
}

void evaluate(byte cbody, int year, byte m, byte D, float UT =12, float h=0, float longit=LONGITUDE, float lat=LATITUDE)
{
    if (cbody !=3)
    {
        oe2radec(cbody, year, m, D, UT);
        radec2azimuth(UT, longit, lat);
        radec2riseset(h, longit, lat);
        Rise += TZ + DST(year,m,D);
        Set += TZ + DST(year,m,D);
        Rise += (Rise<0)?24:0;//TODO popraw, bo pokazuje z dnia wczesniej/pozniej (kilka minut roznicy)
        Set -= (Set>24)?24:0;
    }
    else     //moon
    {
        float LT = TZ + DST(year,m,D);//w dniu zmiany czasu od polnocy do godziny zmiany czasu zle wyliczam DST
        float moonRise1,moonSet1, moonRise2, moonSet2;
        oe2radec(3, year, m, D, 24);
        radec2azimuth(24, longit, lat);
        radec2riseset(h, longit, lat);
        moonRise2 = Rise + LT;
        moonSet2 = Set + LT;
        oe2radec(3, year, m, D, 12);
        radec2azimuth(12, longit, lat);
        radec2riseset(h, longit, lat);
        moonRise1 = Rise + LT;
        moonSet1 = Set + LT;
        oe2radec(3, year, m, D, 0);
        radec2azimuth(0, longit, lat);
        radec2riseset(h, longit, lat);
        Rise = Rise + LT;
        Set = Set + LT;
        if ((moonSet2<24) && (Set>24)) Set-=24;
        if ((moonSet2<24) && (moonSet1>24)) moonSet1-=24;
        if ((Set<24) && (moonSet2>24)) Set = 33.56 ; //w tym dniu nie ma zachodu Ksiezyca
        else
        {
            if(moonSet2>24)
            {
                Set-=24;
                moonSet1-=24;
                moonSet2-=24;
            }
            /* zastosuje aproksymacje liniowa dla pary punktow
            gdy zachod jest przed 12 to
            moonSet(0)=y1=a*x1+b; moonSet(12)=y2=a*x2+b - wyliczam wsp. a i b
            nastepnie kiedy prosta y=a*x+b przecina siê z prost¹ y=x
            czyli x=b/(1-a), co po przeksztalceniach daje
            x=(y1*x2-y2*x1)/(x2-x1-y2+y1)
            podobnie gdy zachod jest po 12*/
            if(moonSet1<12)
                Set=12*Set/(12-moonSet1+Set);
            else
                Set=12*(moonSet1*2-moonSet2)/(12-moonSet2+moonSet1);
        }
        if (moonRise1<0) moonRise1+=24;
        if ((Rise<0)&&(moonRise2>0)) Rise=33.56;// w tym dniu nie ma wschodu Ksiezyca
        else
        {
            if(Rise<0) Rise+=24;
            if(moonRise2<0) moonRise2+=24;
            if(moonRise1<12)
                Rise=12*Rise/(12-moonRise1+Rise);
            else
                Rise=12*(moonRise1*2-moonRise2)/(12-moonRise2+moonRise1);
        }
        oe2radec(3, year, m, D, UT);
        radec2azimuth(UT, longit, lat);
    }

}

#ifdef ARDU

void setup()
{
    Serial.begin(9600);

}
#define COUT(h) Serial.print(" "); Serial.print(#h); Serial.print(":"); Serial.print((int)floor(h)); \
Serial.print("h"); Serial.print((int)((h-floor(h))*60)) ; Serial.print("m");
#define COUT2(h) Serial.print(" "); Serial.print(#h); Serial.print(":"); Serial.print(h);
void loop()
{
    int y=2018,m=5;
    //cbody 0 - Sun, 1 - Mercury, 2 - Venus, 3 - Moon, 4 - Mars, 5 - Jupiter, 6 - Saturn, 7 - Uran, 8 - Neptun
    byte cbody = 5;
    for(int i=1; i<=31; i++)
    {
        Serial.print("<");
        Serial.print(((long)y*100+m)*100+i);
        Serial.print(">");
        evaluate(cbody, y, m, i, 17, 0);
        COUT2(cbody)
        COUT2(RA)
        COUT2(Dec)
        COUT(Rise)
        COUT(Set)
        COUT2(Azimuth)
        COUT2(Altitude)
        Serial.println();
    }
    delay(10000);
}
#else

#define COUT(h) cout <<" " << #h << ":"<<floor(h) << "h" << (int)((h-floor(h))*60) << "m";
#define COUT2(h) cout <<" " << #h << ":"<<(float)h;
int main()
{
    int y=2018,m=5;
    byte cbody =3;
    for(int i=1; i<=31; i++)
    {
        cout << "<"<< (y*100+m)*100+i <<">";
        evaluate(cbody, y, m, i, 12, 0);
        COUT2(cbody)
        COUT2(RA)
        COUT2(Dec)
        COUT(Rise)
        COUT(Set)
        COUT2(Azimuth)
        COUT2(Altitude)
        cout<<endl;
    }
    return 0;
}
#endif // ARDU

