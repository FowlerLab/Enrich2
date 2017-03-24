#  Copyright 2016 Alan F Rubin
#
#  This file is part of Enrich2.
#
#  Enrich2 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Enrich2 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Enrich2.  If not, see <http://www.gnu.org/licenses/>.

def pretty_class_str(obj):
    cls_name = str(obj.__class__).split('.')[-1][:-2]
    string = "{}\n{}\n".format(cls_name, '-'*len(cls_name))
    for (attr, value) in obj.__dict__.items():
        if not(attr.startswith("__") and attr.endswith("__")):
            string += "{}={}\n".format(attr, value)
    return string


    
