#ifndef _CONSOLE_BUF_HPP_
#define _CONSOLE_BUF_HPP_

#include <string>
#include <streambuf>
#include <iomanip>
#include <cstdio>
#include <ctime>

class console_buf : public std::streambuf
{
	protected:
		bool _print_new_line;
		std::time_t _start;
		std::time_t _t;
		struct tm* _now;
		int _current_minute;
		int _time_difference;

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

					printf("%3d:%c%d:%c%d:  ", _time_difference / 3600, (elapsed_minute < 10 ? '0' : 0), elapsed_minute, (elapsed_second < 10 ? '0' : 0), elapsed_second);

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
		console_buf() :
			_print_new_line(true),
			_start(time(0)),
			_t(_start),
			_now(localtime(&_t)),
			_current_minute(-1) {}

};

// (c) 2014 Jerry Coffin
// Graciously borrowed from http://codereview.stackexchange.com/questions/54845/filtering-streambuf
class widthbuf : public std::streambuf
{
	public:
		widthbuf(int w, std::streambuf* s) :
			indent_width(0),
			def_width(w),
			width(w),
			sbuf(s),
			count(0) {}

		~widthbuf() { overflow('\n'); }

		void set_indent(int w)
		{
			if (w == 0)
			{
				prefix.clear();
				indent_width = 0;
				width = def_width;
			}
			else
			{
				indent_width += w;
				prefix = std::string(indent_width, space);
				width -= w;
			}
		}

	private:
		// This is basically a line-buffering stream buffer.
		// The algorithm is:
		// - Explicit end of line ("\r" or "\n"): we flush our buffer
		//   to the underlying stream's buffer, and set our record of
		//   the line length to 0.
		// - An "alert" character: sent to the underlying stream
		//   without recording its length, since it doesn't normally
		//   affect the a appearance of the output.
		// - tab: treated as occupying `tab_width` characters, but is
		//   passed through undisturbed (but if we wanted to expand it
		//   to `tab_width` spaces, that would be pretty easy to do so
		//   you could adjust the tab-width if you wanted.
		// - Everything else: really basic buffering with word wrapping.
		//   We try to add the character to the buffer, and if it exceeds
		//   our line width, we search for the last space/tab in the
		//   buffer and break the line there. If there is no space/tab,
		//   we break the line at the limit.
		int_type overflow(int_type c)
		{
			if (traits_type::eq_int_type(traits_type::eof(), c))
				return traits_type::not_eof(c);

			switch (c)
			{
				case '\n':
				case '\r':
				{
					buffer += c;
					count = 0;
					sbuf->sputn(prefix.c_str(), indent_width);
					int_type rc = sbuf->sputn(buffer.c_str(), buffer.size());
					buffer.clear();
					return rc;
				}
				case '\a':
					return sbuf->sputc(c);
				case '\t':
					buffer += c;
					count += tab_width - count % tab_width;
					return c;
				default:
					if (count >= width)
					{
						size_t wpos = buffer.find_last_of(" \t");
						sbuf->sputn(prefix.c_str(), indent_width);

						if (wpos != std::string::npos)
						{
							sbuf->sputn(buffer.c_str(), wpos);
							count = buffer.size() - wpos - 1;
							buffer = std::string(buffer, wpos + 1);
						}
						else
						{
							sbuf->sputn(buffer.c_str(), buffer.size());
							buffer.clear();
							count = 0;
						}

						sbuf->sputc('\n');
					}

					buffer += c;
					++count;
					return c;
			}
		}

		size_t indent_width;
		size_t width, def_width;
		size_t count;
		size_t tab_count;
		static const int tab_width = 8;
		char_type space = static_cast<char_type>(' ');

		std::string prefix;
		std::string buffer;
		std::streambuf* sbuf;
};

class widthstream : public std::ostream
{
	widthbuf buf;

	public:
		widthstream(size_t width, std::ostream& os) :
			buf(width, os.rdbuf()),
			std::ostream(&buf) {}

		widthstream& indent(int w)
		{
			buf.set_indent(w);
			return *this;
		}

};

#endif	// _CONSOLE_BUF_HPP_