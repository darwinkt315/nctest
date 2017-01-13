#include <iostream>
#include <ncurses.h>
#include <vector>
#include <array>
#include <chrono>
#include <thread>
#include <random>
#include <cassert>
#include <math.h>

using namespace std;

//======================================================================================================

template <class T>
class Rand {
private:
    static default_random_engine gen;
public:
    T operator()(T min, T max) {
        assert(false);
    }
};
template <class T>
default_random_engine Rand<T>::gen(chrono::system_clock::now().time_since_epoch().count());

template<> int Rand<int>::operator()(int min, int max) {
    assert(min<max);
    static uniform_int_distribution<int> uid;
    return uid(gen) % (max-min+1) + min;
}

template<> uint Rand<uint>::operator()(uint min, uint max) {
    assert(min<max);
    static uniform_int_distribution<uint> uid;
    return uid(gen) % (max-min+1) + min;
}

template<> char Rand<char>::operator()(char min, char max) {
    assert(min<max);
    static uniform_int_distribution<char> uid;
    return uid(gen) % (max-min+1) + min;
}

template<> double Rand<double>::operator()(double min, double max) {
    assert(min<max);
    static uniform_real_distribution<double> uid;
    return uid(gen) * (max-min) + min;
}

template<> float Rand<float>::operator()(float min, float max) {
    assert(min<max);
    static uniform_real_distribution<float> uid;
    return uid(gen) * (max-min) + min;
}

//======================================================================================================

template <class T>
struct Point {
    T x;
    T y;
};

typedef Point<int>      IPoint;
typedef Point<uint>     UPoint;
typedef Point<char>     CPoint;
typedef Point<double>   DPoint;
typedef Point<float>    FPoint;

//======================================================================================================

template <class T>
struct Range {
private:
    static Rand<T> rnd;
public:
    T min;
    T max;
    inline bool ContainsStrict(const T val) const { return min<val && val<max; }
    inline bool Contains(const T val) const { return min<=val && val<=max; }
    inline T Any() const { return rnd(min, max); }
};
template <class T>
Rand<T> Range<T>::rnd;

typedef Range<int>      IRange;
typedef Range<uint>     URange;
typedef Range<char>     CRange;
typedef Range<double>   DRange;
typedef Range<float>    FRange;

//======================================================================================================

typedef  double (*SequenceFunc)(double);

// equation - y = f(x), where x[0..1], y[0..1]
// size => result.size()
// minY => result.min()
// maxY => result.max()
vector<int> GetSequence(uint size, int minY, int maxY, SequenceFunc sf) {
    vector<int> result(size);
    double step = 1.0/(size-1);
    double x = 0;
    int scale = maxY-minY;
    for(uint i = 0; i<size ; ++i, x+=step)
        result[i] = round(scale * sf(x)) + minY;
    return result;
}

double Sin(double x) { return (sin(x*2*M_PI) + 1)/2; }
double Linear(double x) { return x; }
double Cavity(double x) {
    if(x<0.5) return 1-x*2;
    else if(x>0.5) return (x-0.5)*2;
    return 0;
}
double Peak(double x) {
    if(x<0.5) return x*2;
    else if(x>0.5) return (1-x)*2;
    return 1;
}

double BellShaped(double x) {

    const static double ln=1.0/3.0;
    const static double x1 = 0.5 - ln/2;
    const static double x2 = 0.5 + ln/2;

    if(x<x1) return x/x1;
    else if(x1<=x && x<=x2) return 1;
    else return (1-x)/x1;
}

//======================================================================================================

class Timer {
protected:
    uint tickLimit;
    uint tickCurrent;
    bool finishFlag;

public:
    Timer(uint _tickLimit = 0, uint tickOffset=0) :
        tickLimit(_tickLimit), tickCurrent(tickOffset), finishFlag(false) { }
    bool Tick() {
        if(finishFlag) return true;
        if(++tickCurrent > tickLimit) finishFlag = true;
        return finishFlag;
    }
    void Reset(uint offset = 0) { tickCurrent = offset<tickLimit ? offset : 0; finishFlag = false; }
    Timer& operator=(uint value) {
        Reset();
        tickLimit = value;
        return *this;
    }
    uint TickLimit() const      { return tickLimit; }
    uint TickCurrent() const    { return tickCurrent; }
    bool IsFinish() const       { return finishFlag; }
};

//======================================================================================================

class CycleTimer : public Timer {
protected:
    uint cycleLimit;
    uint cycleCurrent;
public:
    CycleTimer(uint _tickLimit = 0, uint tickOffset = 0, uint _cycleLimit = 0) :
        Timer(_tickLimit, tickOffset), cycleLimit(_cycleLimit), cycleCurrent(0) { }

    bool Tick() {
        if(finishFlag) return true;
        if(++tickCurrent>tickLimit)
        {
            tickCurrent = 0;
            if(cycleLimit)
                if(++cycleCurrent>=cycleLimit)
                    finishFlag = true;
        }
        return finishFlag;
    }

    void Reset(uint offset=0) {
        Timer::Reset(offset);
        cycleCurrent = 0;
    }

    CycleTimer& operator=(uint value) {
        Reset();
        tickLimit = value;
        return *this;
    }

    uint CycleCurrent() const   { return cycleCurrent; }
    uint CycleLimit() const     { return cycleLimit; }
};

//======================================================================================================

class Gradient {
//-----------------------fields-------------------------------
private:
    enum {
        MAX_NC_BR = 1000,
        START_NUM = 20,
        STEPS = 21,
        BG_COLOR = COLOR_BLACK,
        BR_STEP = MAX_NC_BR/(STEPS-1)
    };

public:
    enum class IS {
        GREEN2BLACK,
        GREEN2WHITE,
        WHITE2BLACK
    };

    enum {
        MAX_BRIGHT = STEPS-1
    };

    static short green2black[STEPS];
    static short green2white[STEPS];
    static short white2black[STEPS];

    short begin;
    short level;

//-----------------------methods------------------------------
private:
    inline void CheckBounds() { if(level<0) level=0; else if(level>MAX_BRIGHT) level=MAX_BRIGHT; }

public:
    Gradient(IS val = IS::WHITE2BLACK) : begin{}, level{} {
        switch (val) {
        case IS::GREEN2BLACK:
             begin = green2black[0];
            break;
        case IS::GREEN2WHITE:
            begin = green2white[0];
            break;
        case IS::WHITE2BLACK:
            begin = white2black[0];
            break;
        default:

            break;
        }
    }

    inline void Level(short lvl) { level = lvl; CheckBounds(); }
    inline short Level() const { return level; }

    inline void operator+=(short val) { level+=val; CheckBounds(); }
    inline void operator-=(short val) { level-=val; CheckBounds(); }

    inline Gradient& operator++() { ++level; CheckBounds(); return *this; }
    inline Gradient& operator--() { --level; CheckBounds(); return *this; }

    inline short clr() const { return begin+level; }

    static void Init() {
        //green2black
        for(uint i=START_NUM, j=MAX_NC_BR, k=0; i<START_NUM+STEPS; ++i, j-=BR_STEP, ++k) {
            init_color(i, 0, j, 0);
            init_pair(i, i, BG_COLOR);
            green2black[k] = i;
        }
        //green2white
        for(uint i=START_NUM+STEPS, j=0, k=0; i<START_NUM+STEPS*2; ++i, j+=BR_STEP, ++k) {
            init_color(i, j, 1000, j);
            init_pair(i, i, BG_COLOR);
            green2white[k] = i;
        }
        //white2black
        for(uint i=START_NUM+STEPS*2, j=MAX_NC_BR, k=0; i<START_NUM+STEPS*3; ++i, j-=BR_STEP, ++k) {
            init_color(i, j, j, j);
            init_pair(i, i, BG_COLOR);
            white2black[k] = i;
        }
    }

    static void Test() {
        auto coloredPrint = [](const char* str ,const Gradient& g) {
            attron(COLOR_PAIR(g.clr()));
            printw(str);
            attroff(COLOR_PAIR(g.clr()));
        };

        Gradient g1 = Gradient::IS::GREEN2BLACK;
        Gradient g2 = Gradient::IS::GREEN2WHITE;
        Gradient g3 = Gradient::IS::WHITE2BLACK;

        for(short i=0; i<=STEPS; ++i, ++g1, ++g2, ++g3){
            coloredPrint("###  ", g1);
            coloredPrint("###  ", g2);
            coloredPrint("###\n", g3);
            //printw("%d  %d  %d\n", g1.clr(), g2.clr(), g3.clr());
        }

    }
};
short Gradient::green2black[Gradient::STEPS];
short Gradient::green2white[Gradient::STEPS];
short Gradient::white2black[Gradient::STEPS];

//======================================================================================================
vector<int> GetGradientSequence(uint size, SequenceFunc sf) {
    return GetSequence(size, 0, Gradient::MAX_BRIGHT, sf);
}

class TimedSequence : public CycleTimer {
protected:
    const vector<int>* sequence;
public:
    TimedSequence(const vector<int>* _sequence = nullptr, uint _tickOffset = 0, uint _cycleLimit = 0) :
        CycleTimer((_sequence ? _sequence->size()-1 : 0), _tickOffset, _cycleLimit ), sequence(_sequence) { }

    TimedSequence& operator=(uint value) { return *this; }

    int ElemCurrent() const { return sequence->operator [](tickCurrent); }
    const vector<int>& Sequence() const { return *sequence; }
    void Sequence(const vector<int>& _sequance) { sequence = &_sequance; tickLimit = sequence->size()-1; }

    static void Test() {
        auto sq = GetGradientSequence(21, Cavity);
        TimedSequence ts(&sq,0,1);
        for(uint i=0; i<25; ++i) {
            cout << ts.ElemCurrent() << endl;
            ts.Tick();
        }
    }
};

//======================================================================================================

struct Symbol {
    char ch;
    short clr;
};

//======================================================================================================

namespace Colors {
    enum {
        MAX_BR = 1000,
        START_NUM = 20,
        STEPS = 20
    };

    short green2black[STEPS];
    short green2white[STEPS];
    short white2black[STEPS];

    short bgClr = COLOR_BLACK;

    uint brStep = MAX_BR/STEPS;

    void Init() {
        //green2black
        for(uint i=START_NUM, j=MAX_BR, k=0; i<START_NUM+STEPS; ++i, j-=brStep, ++k) {
            init_color(i, 0, j, 0);
            init_pair(i, i, bgClr);
            green2black[k] = i;
        }
        //green2white
        for(uint i=START_NUM+STEPS, j=0, k=0; i<START_NUM+STEPS*2; ++i, j+=brStep, ++k) {
            init_color(i, j, 1000, j);
            init_pair(i, i, bgClr);
            green2white[k] = i;
        }
        //white2black
        for(uint i=START_NUM+STEPS*2, j=MAX_BR, k=0; i<START_NUM+STEPS*3; ++i, j-=brStep, ++k) {
            init_color(i, j, j, j);
            init_pair(i, i, bgClr);
            white2black[k] = i;
        }
    }

    void Test() {
        auto coloredPrint = [](const char* str, short clr) {
            attron(COLOR_PAIR(clr));
            printw(str);
            attroff(COLOR_PAIR(clr));
        };

        for(uint i=0; i<STEPS; ++i){
            coloredPrint("###  ", green2black[i]);
            coloredPrint("###  ", green2white[i]);
            coloredPrint("###\n", white2black[i]);
        }
    }
}

//======================================================================================================

struct Area {
    vector<UPoint> points;
    vector<char> chars;

    Area() {}
    Area(vector<UPoint> _points, vector<char> _chars) :
        points(_points), chars(_chars) {}
};

//======================================================================================================

typedef vector<UPoint> Track;
class TrackGen {
private:
    static const DRange axisRange;
    static const DRange ksiRange;
    static const URange angleRange;

    UPoint dims;
    DPoint step;
    DPoint pt;

public:
    TrackGen(const UPoint& _dims) : pt{} {
        assert(dims.x && dims.y);
        dims = _dims;
        step = {1.0/dims.x, 1.0/dims.y};
    }

    Track GetLinear(double fillCfnt = 0.8) {
        assert(0<fillCfnt && fillCfnt<1);

        Track result;
        DRange aRange = {0.5 - fillCfnt/2, 0.5 + fillCfnt/2};
        DPoint a = {aRange.Any(), aRange.Any()};
        uint alpha = angleRange.Any();
        DPoint n = {cos(alpha/180.0*M_PI), sin(alpha/180.0*M_PI)};
        pt={0,0};

        if(ksiRange.Contains(alpha)) {
            for(uint i=0; i<dims.y; ++i, pt.y+=step.y) {
                pt.x = n.x/n.y*(pt.y-a.y)+a.x;
                if(axisRange.Contains(pt.x)) result.push_back({uint(pt.x*(dims.x-1)), i});
            }
        }
        else {
            for(uint i=0; i<dims.x; ++i, pt.x+=step.x) {
                pt.y = (n.y/n.x)*(pt.x-a.x)+a.y;
                if(axisRange.Contains(pt.y)) result.push_back({i, uint(pt.y*(dims.y-1))});
            }
        }
        return result;
    }

    Track GetVertical(uint pos) {
        assert(pos<dims.x-1);
        Track result(dims.y);
        for(uint i=0; i<dims.y; ++i)
            result[i] = {pos, i};
        return result;
    }

    Track GetCircle(UPoint center, uint radius) {
        assert(radius > 0);

    }

    static void Test() {
        initscr();
        start_color();
        curs_set(0);
        Gradient::Init();
        UPoint dims;
        getmaxyx(stdscr, dims.y, dims.x);
        vector<UPoint> res;
        TrackGen tg(dims);

        for(uint i=0; i<100; ++i) {
            res = tg.GetLinear();
            for(uint j=0; j<res.size(); ++j)
                mvaddch(res[j].y, res[j].x, '+');

            refresh();
            getch();
            clear();
        }
        refresh();
        getch();
        endwin();
    }
};
const DRange TrackGen::axisRange = {0,1};
const DRange TrackGen::ksiRange = {45, 135};
const URange TrackGen::angleRange = {40, 140};

//======================================================================================================

class Chain {
//-------------------------fields------------------------------
private:
    static const URange chainLengthRange;
    static const CRange charsRange;
    static const URange initDelayRange;
    static const URange moveDelayRange;
    static const uint   TAIL_SIZE = 10;

    string          symbols;
    vector<uint>    colors;
    int             position;
    Track           track;
    bool            finishFlag;
    bool            movedFlag;

    Timer           mainTimer;
    Timer           initDelay;
    Timer           moveDelay;
public:

//-------------------------methods-----------------------------
private:
public:
    Chain(Track&& _track) {
        track = _track;
        Init();
    }

    void Init() {
        finishFlag = false;
        position = -1;
        symbols.resize(chainLengthRange.Any());
        for(uint i=0; i<symbols.size(); ++i) symbols[i] = charsRange.Any();

        colors.resize(symbols.size());
        colors[0] = Gradient::white2black[0];
        colors[1] = Gradient::green2white[Gradient::MAX_BRIGHT/2];
        for(uint i=2; i<colors.size()-TAIL_SIZE; ++i)
            colors[i] = Gradient::green2black[0];
        for(uint i=colors.size()-TAIL_SIZE, j=0; i<colors.size(); ++i, j+=Gradient::MAX_BRIGHT/TAIL_SIZE)
            colors[i] = Gradient::green2black[j];

        mainTimer = symbols.size()+track.size()+1;
        moveDelay = moveDelayRange.Any();
        initDelay = initDelayRange.Any();
    }

    void Init(Track&& _track) {
        track = _track;
        Init();
    }

    void SetTrack(Track&& _track) { track = _track; }
    const Track& GetTrack() const { return track; }

    bool Tick() {
        if(finishFlag) return true;
        if(!initDelay.Tick()) return false;

        if(!moveDelay.Tick()) return false;
        else {
            moveDelay.Reset();
            ++position;
            movedFlag = true;
            if(mainTimer.Tick()) finishFlag=true;
        }
        return finishFlag;
    }

    void Display() {

        if(finishFlag || !movedFlag) return;

        int chainStart = max(0, int(position-track.size())+1);
        int chainEnd = min(position, int(symbols.size())-1);
        int trackStart = max(0, int(position-symbols.size())+1);
        int trackEnd = min(position, int(track.size())-1);

        if(trackStart>0) mvaddch(track[trackStart-1].y, track[trackStart-1].x, ' ');

        for(; chainStart<=chainEnd && trackEnd>=trackStart; ++chainStart, --trackEnd) {
            attron(COLOR_PAIR(colors[chainStart]));
            mvaddch(track[trackEnd].y, track[trackEnd].x, symbols[chainStart]);
            attroff(COLOR_PAIR(colors[chainStart]));
        }
        movedFlag = false;
    }

    static void VerticalTest() {
        initscr();
        start_color();
        curs_set(0);
        nodelay(stdscr, true);
        Gradient::Init();
        UPoint dims;
        getmaxyx(stdscr, dims.y, dims.x);

        TrackGen tg(dims);

        vector<Chain*> mx(dims.x/2);
        for(uint i=0; i<dims.x/2; ++i)
            mx[i] = new Chain(tg.GetVertical(i*2));

        while(true) {

            for(uint i=0; i<mx.size(); ++i) {
                if(mx[i]->Tick()) mx[i]->Init();
                mx[i]->Display();

            }
            refresh();
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
            if(getch() != ERR) break;
        }
        for(uint i=0; i<mx.size(); ++i) delete mx[i];
        endwin();
    }

    static void LinearTest() {
        initscr();
        start_color();
        curs_set(0);
        nodelay(stdscr, true);
        Gradient::Init();
        UPoint dims;
        getmaxyx(stdscr, dims.y, dims.x);

        TrackGen tg(dims);

        static const uint value = 100;

        vector<Chain*> mx(value);
        for(uint i=0; i<mx.size(); ++i)
            mx[i] = new Chain(tg.GetLinear());

        while(true) {

            for(uint i=0; i<mx.size(); ++i) {
                if(mx[i]->Tick()) mx[i]->Init(tg.GetLinear());
                mx[i]->Display();

            }
            refresh();
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
            if(getch() != ERR) break;
        }
        for(uint i=0; i<mx.size(); ++i) delete mx[i];
        endwin();
    }
};
const URange Chain::chainLengthRange = {20, 40};
const CRange Chain::charsRange = {33,126};
const URange Chain::initDelayRange = {20,70};
const URange Chain::moveDelayRange = {20,70};

//======================================================================================================

/*
class Chain : RandomGenerator {
//-------------------------fields------------------------------
private:
    enum { MIN_LENGTH = 70, MAX_LENGTH = 120 }; //PERCENT
    enum { MIN_CHAR = 33, MAX_CHAR = 126 };
    enum { MIN_INIT_DELAY = 20, MAX_INIT_DELAY = 70 };
    enum { MIN_MOVE_DELAY = 20, MAX_MOVE_DELAY = 70 };
    enum { MIN_MUTATE_DELAY = 100, MAX_MUTATE_DELAY = 200 };
    enum STATE {INIT_DELAY, MOVE, MOVE_DELAY };
    enum { MUTATE_CFNT = 10, TAIL_SIZE = 10 };

    vector<Symbol> field;
    UPoint position;
    UPoint dims;
    STATE state;

    Timer initDelay;
    Timer moveDelay;
    Timer mutateDelay;
    Timer gradientHeadDelay;

public:

//------------------------methods------------------------------
private:
    void Init() {
        uint size = GetRandomUnsigned(dims.y*MIN_LENGTH/100, dims.y*MAX_LENGTH/100);
        field.resize(size);
        for(uint i=0; i<size; ++i)
            field[i].ch = char(GetRandomUnsigned(MIN_CHAR, MAX_CHAR));

        position.y = -1;
        initDelay = GetRandomUnsigned(MIN_INIT_DELAY, MAX_INIT_DELAY)/10*10;
        moveDelay = GetRandomUnsigned(MIN_MOVE_DELAY, MAX_MOVE_DELAY)/10*10;
        mutateDelay = GetRandomUnsigned(MIN_MUTATE_DELAY, MAX_MUTATE_DELAY)/10*10;
        gradientHeadDelay = moveDelay.TickLimit()/10;

        field[0].clr = Colors::white2black[Colors::STEPS-1];
        field[1].clr = Colors::green2white[Colors::STEPS-1];
        field[2].clr = Colors::green2white[10];


        for(uint i=3; i<size-TAIL_SIZE; ++i) field[i].clr = Colors::green2black[0];
        for(uint i=1; i<=TAIL_SIZE; ++i) field[size-i].clr =
                Colors::green2black[Colors::STEPS/TAIL_SIZE*(TAIL_SIZE-i)];
    }

    void Mutate() {
        if(state == STATE::INIT_DELAY) return;
        if(mutateDelay.Tick()) return;
        mutateDelay.Reset();

        uint size = MUTATE_CFNT*field.size();
        uint idx = 0;
        for(uint i=0; i<size; ++i) {
            idx = GetRandomUnsigned(0, field.size()-1);
            field[idx].ch = char(GetRandomUnsigned(MIN_CHAR, MAX_CHAR));
        }
    }


public:
    Chain(UPoint _position, UPoint _dims) :
        position(_position), dims(_dims)
    {
        Init();
        state = STATE::INIT_DELAY;
    }

    const UPoint& Position() const { return position; }
    const vector<Symbol>& Field() const { return field; }

    void Tick(){
        //Mutate();
        switch (state) {
        case STATE::INIT_DELAY:
            if (initDelay.Tick()) break;
            ++position.y;
            state = STATE::MOVE_DELAY;
            break;

        case STATE::MOVE: {
            ++position.y;

            field[0].clr = Colors::white2black[Colors::STEPS-1];
            field[1].clr = Colors::green2white[Colors::STEPS-1];
            field[2].clr = Colors::green2white[10];

            if(position.y < dims.y+field.size())
                state = STATE::MOVE_DELAY;
            else {
                Init();
                state = STATE::INIT_DELAY;
            }
            break;
        }

        case STATE::MOVE_DELAY: {

            if(!gradientHeadDelay.Tick()) {
                gradientHeadDelay.Reset();
                field[0].clr = field[0].clr-2 > 0 ? field[0].clr-2 : 0;
                field[1].clr -= 1;
                field[2].clr -= 1;
            }

            if(moveDelay.Tick()) break;
            moveDelay.Reset();
            state = STATE::MOVE;
            break;
        }

        default:
            break;
        }
    }

    void Print() {
        for(uint i=0; i<field.size(); ++i) {
            attron(COLOR_PAIR(field[i].clr));
            mvaddch(position.y-i-1, position.x, ' ');
            mvaddch(position.y-i, position.x, field[i].ch);
            attroff(COLOR_PAIR(field[i].clr));
        }
    }

    static void Test() {
        initscr();
        start_color();
        curs_set(0);
        nodelay(stdscr, true);
        Colors::Init();
        UPoint dims;
        getmaxyx(stdscr, dims.y, dims.x);
        vector<Chain*> vd(dims.x/2);

        for(uint i=0; i<vd.size(); ++i) vd[i] = new Chain({i*2, 0}, {dims.x, dims.y});

        for(uint i=0; i<10000; ++i)
            for(auto item : vd)
                item->Tick();

        while(true) {
            for(auto item : vd) {
                item->Tick();
                item->Print();
            }
            refresh();
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            if(getch() != ERR) break;
        }
        for(auto item : vd) delete item;
        endwin();
    }
};


//======================================================================================================

class Blink : RandomGenerator {
//----------------------------fields------------------------------
    enum { MIN_BLINK_DELAY = 20, MAX_BLINK_DELAY = 100 };
    enum { GRADIENT_INTERVAL = 30, GRADIENT_INTERVAL_OFFSET = 12 };
    enum { MIDDELE_DELAY = GRADIENT_INTERVAL/5, OUTER_DELAY = GRADIENT_INTERVAL/2 };
    enum { CENTER_OFFSET = GRADIENT_INTERVAL/2-5 };

private:
    static const Area centerArea;
    static const Area middleArea;
    static const Area outerArea;

    static const vector<int> sequence;

    TimedSequence centerGrad;
    TimedSequence middleGrad;
    TimedSequence outerGrad;

    Timer mainTimer;

    UPoint position;
    UPoint dims;

public:
//----------------------------methods-----------------------------
private:

public:
    Blink(UPoint _dims) :
        centerGrad(&sequence,CENTER_OFFSET,1),
        middleGrad(&sequence,0,1),
        outerGrad(&sequence,0,1),
        mainTimer(OUTER_DELAY+ GRADIENT_INTERVAL +
                  GetRandomUnsigned(MIN_BLINK_DELAY, MAX_BLINK_DELAY)),
        position{GetRandomUnsigned(0, _dims.x-5), GetRandomUnsigned(0, _dims.y-3)},
        dims(_dims) { }

    void Reset() {
        position = {GetRandomUnsigned(0, dims.x-5), GetRandomUnsigned(0, dims.y-3)};
        centerGrad.Reset(CENTER_OFFSET);
        middleGrad.Reset();
        outerGrad.Reset();
        mainTimer = OUTER_DELAY+ GRADIENT_INTERVAL +
                GetRandomUnsigned(MIN_BLINK_DELAY, MAX_BLINK_DELAY);
    }

    void Tick() {
        if(mainTimer.Tick()) Reset();
        centerGrad.Tick();
        if(mainTimer.TickCurrent() >= MIDDELE_DELAY)
            middleGrad.Tick();
        if(mainTimer.TickCurrent() >= OUTER_DELAY)
            outerGrad.Tick();
    }

    void Print() {
        for(uint i=0; i<centerArea.points.size(); ++i) {
            attron(COLOR_PAIR(Gradient::white2black[centerGrad.ElemCurrent()]));
            mvaddch(centerArea.points[i].y+position.y, centerArea.points[i].x+position.x, centerArea.chars[i]);
            attroff(COLOR_PAIR(Gradient::white2black[centerGrad.ElemCurrent()]));
        }

        for(uint i=0; i<middleArea.points.size(); ++i) {
            attron(COLOR_PAIR(Gradient::white2black[middleGrad.ElemCurrent()]));
            mvaddch(middleArea.points[i].y+position.y, middleArea.points[i].x+position.x, middleArea.chars[i]);
            attroff(COLOR_PAIR(Gradient::white2black[middleGrad.ElemCurrent()]));
        }

        for(uint i=0; i<outerArea.points.size(); ++i) {
            attron(COLOR_PAIR(Gradient::white2black[outerGrad.ElemCurrent()]));
            mvaddch(outerArea.points[i].y+position.y, outerArea.points[i].x+position.x, outerArea.chars[i]);
            attroff(COLOR_PAIR(Gradient::white2black[outerGrad.ElemCurrent()]));
        }
    }

    static void Test() {
        initscr();
        start_color();
        curs_set(0);
        Gradient::Init();
        Colors::Init();
        UPoint dims;
        getmaxyx(stdscr, dims.y, dims.x);
        nodelay(stdscr, true);

        vector<Chain*> vd(dims.x/2);
        for(uint i=0; i<vd.size(); ++i) vd[i] = new Chain({i*2, 0}, {dims.x, dims.y});

        vector<Blink*> bls(dims.x/10);
        for(auto& item: bls) item = new Blink(dims);

        for(uint i=0; i<10000; ++i) {
            for(auto& item : bls)
                item->Tick();
            for(auto item : vd)
                item->Tick();
        }

        while(true) {
            for(auto item : vd) {
                item->Tick();
                item->Print();
            }

            for(auto& item : bls) {
                item->Tick();
                item->Print();
            }
            refresh();
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
            if(getch() != ERR) break;
        }

        for(auto item : vd) delete item;
        for(auto item:bls) delete item;
        endwin();
    }
};
const Area Blink::centerArea{{{2,1}}, {{'X'}}};
const Area Blink::middleArea{{{2,0},{3,1},{2,2},{1,1}},{{'X'},{'X'},{'X'},{'X'}}};
const Area Blink::outerArea{{{3,0},{4,1},{3,2},{1,2},{0,1},{1,0}},{{'X'},{'X'},{'X'},{'X'},{'X'},{'X'}}};
const vector<int> Blink::sequence{GetGradientSequence(Blink::GRADIENT_INTERVAL, Cavity)};
*/

//======================================================================================================

void Greeting() {
    typedef vector<bool> boolVR;
    typedef vector<boolVR> boolMX;

    typedef vector<char> charVR;
    typedef vector<charVR> charMX;

    uint width = 10;
    uint height = 10;
    uint offset = 2;

    static const boolVR empty =
                    {0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,};

    static const boolVR w =
                     {1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,1,1,0,0,1,1,
                     1,1,0,1,1,1,1,0,1,1,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,1,0,0,0,0,1,1,1,
                     1,1,0,0,0,0,0,0,1,1,};

    static const boolVR e =
                    {1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,1,0,0,0,0,
                     1,1,1,1,1,1,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,};

    static const boolVR l =
                    {1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,};

    static const boolVR c =
                    {0,0,1,1,1,1,1,1,1,0,
                     0,1,1,1,1,1,1,1,1,1,
                     1,1,1,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,0,0,0,0,0,0,0,0,
                     1,1,1,0,0,0,0,0,1,1,
                     0,1,1,1,1,1,1,1,1,1,
                     0,0,1,1,1,1,1,1,1,0,};

    static const boolVR o =
                    {0,0,1,1,1,1,1,1,0,0,
                     0,1,1,1,1,1,1,1,1,0,
                     1,1,1,0,0,0,0,1,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,1,0,0,0,0,1,1,1,
                     0,1,1,1,1,1,1,1,1,0,
                     0,0,1,1,1,1,1,1,0,0,};

    static const boolVR m =
                    {1,1,0,0,0,0,0,0,1,1,
                     1,1,1,0,0,0,0,1,1,1,
                     1,1,1,1,0,0,1,1,1,1,
                     1,1,0,1,1,1,1,0,1,1,
                     1,1,0,0,1,1,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,};

    static const boolVR t =
                    {1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,};

    static const boolVR a =
                    {0,0,0,0,1,1,0,0,0,0,
                     0,0,0,1,1,1,1,0,0,0,
                     0,0,1,1,0,0,1,1,0,0,
                     0,1,1,0,0,0,0,1,1,0,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,};

    static const boolVR r =
                    {1,1,1,1,1,1,1,1,0,0,
                     1,1,1,1,1,1,1,1,1,0,
                     1,1,0,0,0,0,0,1,1,0,
                     1,1,1,1,1,1,1,1,0,0,
                     1,1,1,1,1,1,1,0,0,0,
                     1,1,0,0,1,1,0,0,0,0,
                     1,1,0,0,0,1,1,0,0,0,
                     1,1,0,0,0,0,1,1,0,0,
                     1,1,0,0,0,0,0,1,1,0,
                     1,1,0,0,0,0,0,0,1,1,};

    static const boolVR i =
                    {0,0,1,1,1,1,1,1,0,0,
                     0,0,1,1,1,1,1,1,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,0,0,1,1,0,0,0,0,
                     0,0,1,1,1,1,1,1,0,0,
                     0,0,1,1,1,1,1,1,0,0,};

    static const boolVR x =
                    {1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,
                     0,1,1,0,0,0,0,1,1,0,
                     0,0,1,1,0,0,1,1,0,0,
                     0,0,0,1,1,1,1,0,0,0,
                     0,0,0,1,1,1,1,0,0,0,
                     0,0,1,1,0,0,1,1,0,0,
                     0,1,1,0,0,0,0,1,1,0,
                     1,1,0,0,0,0,0,0,1,1,
                     1,1,0,0,0,0,0,0,1,1,};


    static const URange symbolRange = {'0', '1'};
    vector<int> bellSeq = GetGradientSequence(60, BellShaped);
    const uint delay = 20;

    charMX fl;
    charMX sl;

    auto FillCharMX = [width, height, offset](charMX& dst, const boolMX& mask) {
        for(uint k=0; k<mask.size(); ++k) {
            for(uint i=0; i<width; ++i) {
                charVR column(height);
                for(uint j=0; j<height; ++j)
                    column[j] = mask[k][j*width+i] ? symbolRange.Any() : 0;
                dst.push_back(column);
            }
            if(k != mask.size()-1)
                for(uint i=0; i<offset; ++i)
                    dst.push_back(charVR(height, 0));
        }
    };

    auto DisplayCharMX = [&bellSeq, delay] (const charMX& mx, const UPoint& offset) {
        for(uint i=0; i<mx.size()+bellSeq.size()+1; ++i) {
            for(uint k=max(0,int(bellSeq.size()-1-i)),
                     m=max(0,int(i-bellSeq.size()-1)),
                     l=min(int(bellSeq.size()-k),int(mx.size()-m)); l>0; ++k, ++m, --l) {

                attron(COLOR_PAIR(Gradient::green2black[Gradient::MAX_BRIGHT-bellSeq[k]]));
                for(uint n=0; n<mx[m].size(); ++n)
                    if(mx[m][n]) mvaddch(n+offset.y, m+offset.x, mx[m][n]);
                attroff(COLOR_PAIR(Gradient::green2black[Gradient::MAX_BRIGHT-bellSeq[k]]));
            }

            refresh();
            std::this_thread::sleep_for(std::chrono::milliseconds(delay));
        }
    };

    FillCharMX(fl, {w,e,l,c,o,m,e,empty,t,o});
    FillCharMX(sl, {m,a,t,r,i,x});

    UPoint margin1 = {5,3};
    UPoint margin2 = {margin1.x+(width+offset)*2, margin1.y+height+offset};

    initscr();
    start_color();
    curs_set(0);
    Gradient::Init();

    DisplayCharMX(fl, margin1);
    DisplayCharMX(sl, margin2);

    refresh();
    getch();
    endwin();
}

//======================================================================================================
int main()
{
    //TimedSequence::Test();
    //Blink::Test();
    //Greeting();
    TrackGen::Test();
    //Chain::VerticalTest();
    //Chain::LinearTest();
}





























