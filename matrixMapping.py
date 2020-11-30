

def matrix_mapping(face_1, face_2, min_loc):

    print("Indentifying minimum correlation and mapping correlation back to cube faces...")
    n1 = face_1.shape[1] - 1
    n2 = face_2.shape[1] - 1
    m1 = face_1.shape[0] - 1
    m2 = face_2.shape[0] - 1

    if n1 >= n2:
        if min_loc[1] >= n1:
            matrix_1_xint = (min_loc[1] - n2,n1 + 1)
            matrix_2_xint = (0,n2 - min_loc[1] + n1 +1)
        else:
            if min_loc[1] <= n2:
                matrix_1_xint = (0, min_loc[1] + 1)
                matrix_2_xint = (n2 - min_loc[1], n2 + 1)
            else:
                matrix_1_xint = (min_loc[1] - n2, min_loc[1] + 1)
                matrix_2_xint = (0, n2 + 1)
    else:
        if min_loc[1] >= n1:
            if min_loc[1] <= n2:
                matrix_1_xint = (0, n1 + 1)
                matrix_2_xint = (n2 - min_loc[1], n2 - min_loc[1] + n1 + 1)
            else:
                matrix_1_xint = (min_loc[1] - n2,n1 + 1)
                matrix_2_xint = (0,n2 + n1 - min_loc[1] + 1)
        else:
            matrix_1_xint = (0, n1 + n2 - min_loc[1] + 1)
            matrix_2_xint = (n2 - min_loc[1], n2 + 1)

    if m1 >= m2:
        if min_loc[0] >= m1:
            matrix_1_yint = (min_loc[0] - m2,m1 + 1)
            matrix_2_yint = (0,m2 - min_loc[0] + m1 + 1)
        else:
            if min_loc[0] <= m2:
                matrix_1_yint = (0, min_loc[0] + 1)
                matrix_2_yint = (m2 - min_loc[0], m2 + 1)
            else:
                matrix_1_yint = (min_loc[0] - m2, min_loc[0] + 1)
                matrix_2_yint = (0, m2 + 1)
    else:
        if min_loc[0] >= m1:
            if min_loc[0] <= m2:
                matrix_1_yint = (0, m1 + 1)
                matrix_2_yint = (m2 - min_loc[0], m2 - min_loc[0] + m1 + 1)
            else:
                matrix_1_yint = (min_loc[0] - m2,m1 + 1)
                matrix_2_yint = (0, m2 +  m1 - min_loc[0] + 1)
        else:
            matrix_1_yint = (0, m1 + m2 - min_loc[0] + 1)
            matrix_2_yint = (m2 - min_loc[0],m2 +1)


    matrix_1 = face_1[matrix_1_yint[0]:matrix_1_yint[1], matrix_1_xint[0]:matrix_1_xint[1]]
    matrix_2 = face_2[matrix_2_yint[0]:matrix_2_yint[1], matrix_2_xint[0]:matrix_2_xint[1]]
    
    
    return matrix_1, matrix_2
