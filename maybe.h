#ifndef MAYBE_H
#define MAYBE_H

/**
 * @class MayBe implements value that can be invalid
 */
template<class T> class MayBe
{
public:
    MayBe()
        :
          mIsValid(false)
    {}

    MayBe(const T& val)
        :
          mIsValid(true),
          mVal(val)
    {}

    operator T() const{ return mVal; }
    operator bool() const{ return mIsValid; }

    const T& val() const{ return mVal; }
    void val(const T& val) { mVal = val; }

    bool isValid() const { return mIsValid; }
    void isValid(bool isValid) { mIsValid = isValid; }

private:
    bool mIsValid;
    T mVal;
};

#endif // MAYBE_H
