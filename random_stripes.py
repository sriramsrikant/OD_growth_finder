import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

'''TODO:
create user input
create checkboxes for color selection
add 2nd color option for rows
add defined 2-color fair isle patterns
'''

# input from user
requested_colors = ['red', 'blue', 'black', 'white']
second_color = 'some'
requested_weights = [1, 1, 1, 4]
requested_widths = [2, 4, 6, 8]
total_rows = 150

n = len(requested_colors)
r = len(requested_widths)

colors = list()  # create list of colors with num of entries reflecting weights
for i in range(n):
    add_col = [requested_colors[i]]*requested_weights[i]
    colors.extend(add_col)
m = len(colors)

# initialize variables
stripes = list()
second_stripes = list()
per_color_rows = [0] * n
per_color_sec_rows = [0] * n  # for second colors
num_rows = 0
curr_col = 'blank'
old_col = 'blank'

while num_rows < total_rows:
    while curr_col == old_col:  # don't have two stripes in a row of same color
        i = np.random.randint(m)
        curr_col = colors[i]

        second_col = curr_col
        if second_color is 'all':  # if option for two-color stripes
            while second_col == curr_col:
                second_col = colors[np.random.randint(m)]  # get second col

        elif second_color is 'some':
            j = np.random.random() # determine whether to have a second color
            if j < 0.5:  # can change this threshold to determine likelihood of 2-color row
                while second_col == curr_col:
                    second_col = colors[np.random.randint(m)]
            else:
                second_col = curr_col

    j = np.random.randint(r)  # get random number for stripe width
    stripe_width = requested_widths[j]

    stripes.append((curr_col, stripe_width))  # add tuple of color, width to stripe list

    if second_color:
        second_stripes.append((second_col, stripe_width))

    num_rows = num_rows + stripe_width  # keep track of rows created
    old_col = curr_col  # reassign for next round

clipped_rows =  num_rows - total_rows  # don't have more than requested num of rows
last_stripe = stripes[-1][1] - clipped_rows  # remove extra rows from last stripe
stripes[-1] = (stripes[-1][0], last_stripe)

for col, wid in stripes:  # calculate per_color_rows
    i = requested_colors.index(col)  # get loc of col in requested_colors list
    per_color_rows[i] += wid  # add wid to that loc in per_color_rows

if second_color:
    for col, wid in second_stripes:  # calculate per_color_rows
        i = requested_colors.index(col)  # get loc of col in requested_colors list
        per_color_sec_rows[i] += wid  # add wid to that loc in per_color_rows

my_cmap = matplotlib.colors.ListedColormap(requested_colors, name='req_cols')
col_array = np.zeros(num_rows)
sec_col_array = np.zeros(num_rows)
i = 0
for col, wid in stripes:
    w = 0
    while w < wid:
        c = requested_colors.index(col)
        col_array[i] = c
        i += 1
        w += 1

if second_color:
    i = 0
    for col, wid in second_stripes:
        w = 0
        while w < wid:
            c = requested_colors.index(col)
            sec_col_array[i] = c
            i += 1
            w += 1
    pcol_array = np.column_stack([col_array, sec_col_array])

else:
    pcol_array = np.column_stack([col_array, col_array])

ax = plt.axes()
plt.pcolor(pcol_array, cmap=my_cmap)
ax.set_ylim(0,num_rows)
ax.yaxis.tick_right()
ax.xaxis.set_visible(False)
plt.title('random stripe pattern')
plt.show()

col_string = [col +': '+str(wid) for col, wid in stripes]
print ', '.join(col_string)

col_totals = list()
for i in range(len(per_color_rows)):
    tot_string = requested_colors[i] +': '+ str(per_color_rows[i])
    col_totals.append(tot_string)
print ', '.join(col_totals)