//将四元数转化为旋转矩阵
void q2tr(double q[4],double r[3][3])
{
    double s,x,y,z;
    s = q[0];
    x = q[1];
    y = q[2];
    z = q[3];

    r[0][0] = 1 - 2 * ( pow(y,2) + pow(z,2));
    r[0][1] = 2 * ( x * y - s * z );
    r[0][2] = 2 * ( x * z + s * y );

    r[1][0] = 2 * ( x * y + s * z );
    r[1][1] = 1 - 2 * ( pow(x,2) + pow(z,2));
    r[1][2] = 2 * ( y * z - s * x );

    r[2][0] = 2 * ( x * z - s * y );
    r[2][1] = 2 * ( y * z + s * x );
    r[2][2] = 1 - 2 * ( pow(x,2) + pow(y,2));
}

//将旋转矩阵转化为四元数
void tr2q(double q[4],double t[3][3])
{
    double w2,x2,y2,z2;
    double w,x,y,z;

    w2 = (t[0][0] + t[1][1] + t[2][2] + 1) / 4.0;

    if(w2 > 0.00001)
    {
        w = sqrt(w2);
        x = (t[1][2] - t[2][1]) /(4.0*w);
        y = (t[2][0] - t[0][3]) /(4.0*w);
        z = (t[0][1] - t[1][0]) /(4.0*w);
    }
    else
    {
        w = 0;
        x2 = (t[0][0] - t[1][1] - t[2][2] + 1) / 4.0;
        if(x2 > 0.00001)
        {
            x = sqrt(x2);
            y = (t[0][1] + t[1][0])/(4.0 * x );
            z = (t[0][2] + t[2][0])/(4.0 * x );
        }
        else
        {
            x = 0;
            y2 = (t[1][1] - t[0][0] - t[2][2] + 1) / 4.0;
            if(y2 > 0.00001)
            {
                y = sqrt(y2);
                z = (t[1][2] + t[2][1]) /(4.0*y);
            }
            else
            {
                y = 0;
                z = 1;
            }
        }
    }
    q[0] = w;
    q[1] = x;
    q[2] = y;
    q[3] = z;
}

//k1 = sin( ( 1 - t ) * theta ) / sin(theta);
//k2 = sin( t * theta ) / sin(theta);
//球面插值公式q = q1 * k1 + q2 * k2;
void Slerp(double q1[4],double q2[4],double q[4],double t)
{
    double theta;
    double k1,k2;
    double cosTheta = 0;
    int i;

    for( i = 0 ; i < 4 ; i++ )
        cosTheta += q1[i] * q2[i];

    if( cosTheta < 0 )
    {
        cosTheta = -cosTheta;
        for( i = 0 ; i < 4; i++ )
            q1[i] = -q1[i];
    }

    theta = acos(cosTheta);

    //0.00001为系统的计算精度
    if((theta <= 0.00001) & (theta >= -0.00001))
    {
        for( i = 0 ; i < 4 ; i++ )
            q[i] = q1[i];
    }
    else
    {
        k1 = sin( ( 1 - t ) * theta ) / sin(theta);
        k2 = sin( t * theta ) / sin(theta);
        for( i = 0 ; i < 4 ; i++ )
            q[i] = q1[i] * k1 + q2[i] * k2;
   }
 
}

//StartMat为起点的姿态
//EndMat为终点的姿态
//Middle为插补的中间姿态
//t为当前的距离
void QuaternionInter(double StartMat[3][3],double EndMat[3][3],double MiddleMat[3][3],double t)
{
    double q1[4],q2[4],q[4];

    //将矩阵转化为四元数
    tr2q(q1,StartMat);
    tr2q(q2,EndMat);
 
    //两点的姿态插补
    Slerp(q1,q2,q,t);
    
    //将四元数转化为矩阵           
    q2tr(q,MiddleMat);    
}
