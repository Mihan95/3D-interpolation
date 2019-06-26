#include <QThread>

class Solve_thread : public QThread
{
    Q_OBJECT
public:
    Solve_thread(){}

    void run()
    {
        exec();
    }
};
