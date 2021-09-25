#ifdef GL_ES 
precision mediump float;
uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;
#else 
uniform vec2 u_mouse;
uniform vec2 u_resolution;
uniform vec3 left_top_position;
uniform float osg_FrameTime;
varying vec3 basic_prop;;
#define u_time osg_FrameTime
#endif

#define NUM_PARTICLES 30.0
#define GLOW 0.5

const float PI = 3.1415926;

float circle_shape(vec2 position, float radius)
{
    position -= 0.5;

    float color = smoothstep(radius, radius - 0.1, length(position));

    return color;
}
vec3 Orb(vec2 uv, vec3 color, float radius, float offset)
{
    vec2 position = vec2(sin(offset * (u_time + 30.)),
                        cos(offset * (u_time + 30.)));
    
    position *= sin((u_time) - offset) * cos(offset);
    radius = radius;

    float dist = radius / distance(uv, position);

    return color * pow(dist, 1.0 / GLOW);
}
float set_top(float loop_count, float particle_delta)
{
    return 0.8 + loop_count / particle_delta;
}
float set_center(float particle_num, float width, float left)
{
    return smoothstep(0., NUM_PARTICLES, particle_num) * width - left;
}
float set_flow_time(float t, float loop_count, float particle_time_factor)
{
    return smoothstep(0., t, u_time + 2.0 - loop_count * particle_time_factor);
}
float Particle(vec2 uv, float radius, float particle_num, float loop_count)
{
    float top = set_top(loop_count, 20.);
    float x = set_center(particle_num, 0.6, 0.3);
    float flow_time = set_flow_time(20., loop_count, 0.5);
    float angle = flow_time * PI / 2., tan_angle = 80.0 * PI / 180.;
    float y = smoothstep(0., 1., sin(angle)) * top - loop_count / 20.;

    if(y < 0.) return 0.;

    if(particle_num < NUM_PARTICLES / 2.)
        x -= flow_time * top / tan(tan_angle);
    else 
        x += flow_time * top / tan(tan_angle);

    vec2 position = vec2(x, y);
    float dist = radius / distance(uv, position);

    return pow(dist, 1.0 / GLOW);
}

void circle_main()
{
    vec2 position = gl_FragCoord.xy / u_resolution;

    float color = circle_shape(position, 0.2);

    gl_FragColor = vec4(vec3(color), 1.);
}
void particle_main()
{
    vec2 uv = vec2(gl_FragCoord.xy * 2. - u_resolution.xy) / u_resolution.xy;
    //vec2 uv = vec2(gl_TexCoord[0].x * 2. - 1., gl_TexCoord[0].y * -2. + 1.);
    vec3 pixel = vec3(0,0,0);
    vec3 color = vec3(0,0,0);

    color.r = ((sin(((u_time)) * 0.55) + 1.5) * 0.4);
    color.g = ((sin(((u_time)) * 0.34) + 2.0) * 0.4);
    color.b = ((sin(((u_time)) * 0.31) + 4.5) * 0.3);

    float radius = 0.005;

    for(float j = 0.; j < 10.; ++j)
    {
        for(float i = 0.0; i < NUM_PARTICLES; ++i)
        {
            pixel += color * Particle(uv, radius, i, j);
        }
    }

    //gl_FragColor = vec4(uv, 0.8 + 0.5 * sin(u_time), 1.0);
    gl_FragColor = mix(vec4(uv.xy, 0.8 + 0.5 * sin(u_time), 1.0), vec4(pixel, 1.0), 0.8);
}
void test_main()
{
    //vec2 uv = vec2(gl_FragCoord.xy * 2. - u_resolution.xy) / u_resolution;
    //vec3 uv = vec3(sin(basic_prop));
    
    vec2 uv = vec2(gl_FragCoord.xy / u_resolution.xy);
    gl_FragColor = vec4(uv, 0.0, 1.0);

    //gl_FragColor = vec4(gl_TexCoord[0].x, 1.0 - gl_TexCoord[0].y, 0., 1.0);
}

// comment this string to see each part in full screen
#define BOTH
// uncomment this string to see left part
//#define LEFT

//#define LOW_QUALITY

#define DITHERING

//#define TONEMAPPING

//-------------------
#define pi 3.14159265
#define R(p, a) p=cos(a)*p+sin(a)*vec2(p.y, -p.x)

uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;

// iq's noise
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
	f = f*f*(3.0-2.0*f);
	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
	vec2 rg = texture2D( iChannel0, (uv+ 0.5)/256.0, 0.0 ).yx;
	return 1. - 0.82*mix( rg.x, rg.y, f.z );
}

float fbm( vec3 p )
{
   return noise(p*.06125)*.5 + noise(p*.125)*.25 + noise(p*.25)*.125 + noise(p*.4)*.2;
}

float Sphere( vec3 p, float r )
{
    return length(p)-r;
}

//==============================================================
// otaviogood's noise from https://www.shadertoy.com/view/ld2SzK
//--------------------------------------------------------------
// This spiral noise works by successively adding and rotating sin waves while increasing frequency.
// It should work the same on all computers since it's not based on a hash function like some other noises.
// It can be much faster than other noise functions if you're ok with some repetition.
const float nudge = 4.;	// size of perpendicular vector
float normalizer = 1.0 / sqrt(1.0 + nudge*nudge);	// pythagorean theorem on that perpendicular to maintain scale
float SpiralNoiseC(vec3 p)
{
    float n = -mod(u_time * 0.2,-2.); // noise amount
    float iter = 2.0;
    for (int i = 0; i < 8; i++)
    {
        // add sin and cos scaled inverse with the frequency
        n += -abs(sin(p.y*iter) + cos(p.x*iter)) / iter;	// abs for a ridged look
        // rotate by adding perpendicular and scaling down
        p.xy += vec2(p.y, -p.x) * nudge;
        p.xy *= normalizer;
        // rotate on other axis
        p.xz += vec2(p.z, -p.x) * nudge;
        p.xz *= normalizer;
        // increase the frequency
        iter *= 1.733733;
    }
    return n;
}

float VolumetricExplosion(vec3 p)
{
    float final = Sphere(p,4.);
    #ifdef LOW_QUALITY
    final += noise(p*12.5)*.2;
    #else
    final += fbm(p*50.);
    #endif
    final += SpiralNoiseC(p.zxy*0.4132+333.)*3.0; //1.25;

    return final;
}

float map(vec3 p) 
{
	R(p.xz, 0.008*pi+u_time*0.1);

	float VolExplosion = VolumetricExplosion(p/0.5)*0.5; // scale
    
	return VolExplosion;
}
//--------------------------------------------------------------

// assign color to the media
vec3 computeColor( float density, float radius )
{
	// color based on density alone, gives impression of occlusion within
	// the media
	vec3 result = mix( vec3(1.0,0.9,0.8), vec3(0.4,0.15,0.1), density );
	
	// color added to the media
	vec3 colCenter = 7.*vec3(0.8,1.0,1.0);
	vec3 colEdge = 1.5*vec3(0.48,0.53,0.5);
	result *= mix( colCenter, colEdge, min( (radius+.05)/.9, 1.15 ) );
	
	return result;
}

bool RaySphereIntersect(vec3 org, vec3 dir, out float near, out float far)
{
	float b = dot(dir, org);
	float c = dot(org, org) - 8.;
	float delta = b*b - c;
	if( delta < 0.0) 
		return false;
	float deltasqrt = sqrt(delta);
	near = -b - deltasqrt;
	far = -b + deltasqrt;
	return far > 0.0;
}

// Applies the filmic curve from John Hable's presentation
// More details at : http://filmicgames.com/archives/75
vec3 ToneMapFilmicALU(vec3 _color)
{
	_color = max(vec3(0), _color - vec3(0.004));
	_color = (_color * (6.2*_color + vec3(0.5))) / (_color * (6.2 * _color + vec3(1.7)) + vec3(0.06));
	return _color;
}

void main()
{  
    const float KEY_1 = 49.5/256.0;
	const float KEY_2 = 50.5/256.0;
	const float KEY_3 = 51.5/256.0;
    float key = 0.0;
    key += 0.7*texture2D(iChannel1, vec2(KEY_1,0.25)).x;
    key += 0.7*texture2D(iChannel1, vec2(KEY_2,0.25)).x;
    key += 0.7*texture2D(iChannel1, vec2(KEY_3,0.25)).x;

    //vec2 uv = gl_FragCoord.xy /u_resolution.xy;
    vec2 uv = vec2(gl_TexCoord[0].x, 1.0 - gl_TexCoord[0].y);
    
	// ro: ray origin
	// rd: direction of the ray
	//vec3 rd = normalize(vec3((gl_FragCoord.xy-0.5*u_resolution.xy)/u_resolution.xy, 1.));
    vec3 rd = normalize(vec3(gl_TexCoord[0].x - 0.5, -gl_TexCoord[0].y + 0.5, 1.));
    //vec3 rd = normalize(vec3(()))
	vec3 ro = vec3(0., 0., -6.+key*1.6);
    
	// ld, td: local, total density 
	// w: weighting factor
	float ld=0., td=0., w=0.;

	// t: length of the ray
	// d: distance function
	float d=1., t=0.;
    
    const float h = 0.1;
   
	vec4 sum = vec4(0.0);
   
    float min_dist=0.0, max_dist=0.0;

    if(RaySphereIntersect(ro, rd, min_dist, max_dist))
    {
       
	t = min_dist*step(t,min_dist);
   
	// raymarch loop
    #ifdef LOW_QUALITY
	for (int i=0; i<56; i++)
    #else
    for (int i=0; i<86; i++)
    #endif
	{
	 
		vec3 pos = ro + t*rd;
  
		// Loop break conditions.
	    if(td>0.9 || d<0.12*t || t>10. || sum.a > 0.99 || t>max_dist) break;
        
        // evaluate distance function
        float d = map(pos);
        
        #ifdef BOTH
        /*
        if (uv.x<0.5)
        {
            d = abs(d)+0.07;
        }
        */
        //split screen variant
        //d = uv.x < 0.5 ? abs(d)+0.07 : d;
        
        d = cos(u_time)*uv.x < 0.1 ? abs(d)+0.07 : d;
        #else
        #ifdef LEFT
        d = abs(d)+0.07;
        #endif
		#endif
        
		// change this string to control density 
		d = max(d,0.03);
        
        // point light calculations
        vec3 ldst = vec3(0.0)-pos;
        float lDist = max(length(ldst), 0.001);

        // the color of light 
        vec3 lightColor=vec3(1.0,0.5,0.25);
        
        sum.rgb+=(lightColor/exp(lDist*lDist*lDist*.08)/30.); // bloom
        
		if (d<h) 
		{
			// compute local density 
			ld = h - d;
            
            // compute weighting factor 
			w = (1. - td) * ld;
     
			// accumulate density
			td += w + 1./200.;
		
			vec4 col = vec4( computeColor(td,lDist), td );
            
            // emission
            sum += sum.a * vec4(sum.rgb, 0.0) * 0.2 / lDist;	
            
			// uniform scale density
			col.a *= 0.2;
			// colour by alpha
			col.rgb *= col.a;
			// alpha blend in contribution
			sum = sum + col*(1.0 - sum.a);  
       
		}
      
		td += 1./70.;

        #ifdef DITHERING
        // idea from https://www.shadertoy.com/view/lsj3Dw
        vec2 uvd = uv;
        uvd.y*=120.;
        uvd.x*=280.;
        d=abs(d)*(.8+0.08*texture2D(iChannel2,vec2(uvd.y,-uvd.x+0.5*sin(4.*u_time+uvd.y*4.0))).r);
        #endif 
		
        // trying to optimize step size
        #ifdef LOW_QUALITY
        t += max(d*0.25,0.01);
        #else
        t += max(d * 0.08 * max(min(length(ldst),d),2.0), 0.01);
        #endif
        
        
	}
    
    // simple scattering
    #ifdef LOW_QUALITY    
    sum *= 1. / exp( ld * 0.2 ) * 0.9;
    #else
    sum *= 1. / exp( ld * 0.2 ) * 0.8;
    #endif
        
   	sum = clamp( sum, 0.0, 1.0 );
   
    sum.xyz = sum.xyz*sum.xyz*(3.0-2.0*sum.xyz);
    
	}

    if(length(sum.xy) < 0.1)
        discard;
   
    #ifdef TONEMAPPING
    gl_FragColor = vec4(ToneMapFilmicALU(sum.xyz*2.2),1.0);
	#else
    gl_FragColor = vec4(sum.xyz,1.0);
	#endif
}