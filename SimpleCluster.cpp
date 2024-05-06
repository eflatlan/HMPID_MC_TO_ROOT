#ifndef SimpleCluster_H
#define SimpleCluster_H
class SimpleCluster

{

private:
  double mX, mY, mQ, mSize;

public:
  double x() const { return mX; }

  double y() const { return mY; }

  double q() const { return mQ; }

  double size() const { return mSize; }

  SimpleCluster(double x, double y, double q, double size)

      : mX(x), mY(y), mQ(q), mSize(size) {}

  SimpleCluster() : mX(0), mY(0), mQ(0), mSize(0) {}
};

#endif
