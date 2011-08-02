#ifndef _MagicParser_h
#define _MagicParser_h

class MagicParser
{
public:
  MagicParser() : m_argc(0), m_argv(NULL)
  {
  }

  void SetArgs( int argc, char* argv );

private:
  int   m_argc;
  char* m_argv;
};

#endif
