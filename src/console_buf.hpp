#ifndef _CONSOLE_BUF_HPP_
#define _CONSOLE_BUF_HPP_

#include <streambuf>
#include <cstdio>
#include <ctime>

class console_buf : public std::streambuf
{
protected:
	bool 		_print_new_line;
	std::time_t	_start;
	std::time_t	_t;
	struct tm*	_now;
	int			_current_minute;
	int			_time_difference;
	
	// central formatting function
	virtual int_type overflow (int_type c)
	{
		if (c != EOF)
		{
			// Check whether we are on a new line
			if (_print_new_line)
			{
				// print newline stuff here:
				_t = time(0);
				_now = localtime(&_t);
				
				_time_difference = _t - _start;
				
				int elapsed_minute = (_time_difference / 60) % 60;
				int elapsed_second = _time_difference % 60;
				
				if (_current_minute == _now->tm_min)
				{
					// same minute
					printf("        ");
				}
				else
				{
					// new minute
					char bug[6];
					strftime(bug, 6, "%H:%M", _now);
					printf("%s - ", bug);
				}
				
				printf("%3d:%c%d:%c%d:  ", _time_difference/3600, (elapsed_minute < 10 ? '0' : 0), elapsed_minute, (elapsed_second < 10 ? '0' : 0), elapsed_second);
				
				_current_minute = _now->tm_min;
				_print_new_line = false;
			}
			
			if (c == '\n')
			{
				_print_new_line = true;
			}
			
			// and write the character to the standard output
			if (putchar(c) == EOF)
			{
				return EOF;
			}
		}
		return c;
	}
	
public:
	console_buf() : _print_new_line(true), _start(time(0)), _t(_start), _now(localtime(&_t)), _current_minute(-1) {}
};

#endif // _CONSOLE_BUF_HPP_
