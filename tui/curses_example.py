"""
Example:
https://xerocrypt.wordpress.com/2014/09/03/the-curses-api-and-python/

(C) 2020 Sean Cummins, Drogheda, Ireland
Realeased under the GNU Public License (GPL)
email seancummins16@gmail.com
"""

from os import system
import curses

myscreen = curses.initscr()

myscreen.border(0)
myscreen.addstr(5, 5, "HELLO WORLD")
input = myscreen.getstr(12, 20, 50)
#system("cat /proc/iomem")
myscreen.refresh()
myscreen.getch()

curses.endwin()