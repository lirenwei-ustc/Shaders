#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
//u_resolution 存储了画布的宽高

//画圆
float circle_shape(vec2 position, float radius)
{
    return step(radius, length(position - 0.5));
}
void circle_main()
{
    vec2 position = gl_FragCoord.xy / u_resolution;
    
    vec3 color = vec3(0.);

    float circle = circle_shape(position, 0.2);

    color = vec3(circle);

    gl_FragColor  = vec4(color, 1.);
}

//画长方形
float rect_shape(vec2 position, vec2 scale)
{
    scale = vec2(0.5) - scale * 0.5;

    vec2 shaper = vec2(step(scale.x, position.x), step(scale.y, position.y));

    shaper *= vec2(step(scale.x, 1. - position.x), step(scale.y, 1. - position.y));

    return shaper.x  * shaper.y;
}
void main()
{
    vec2 position = gl_FragCoord.xy / u_resolution;

    vec3 color = vec3(0.);

    float rectangle = rect_shape(position, vec2(0.6, 0.3));
    //width, height

    color = vec3(rectangle);

    gl_FragColor = vec4(color, 1.0);
}