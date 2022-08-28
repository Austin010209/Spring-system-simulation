package comp559.particle;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintStream;


public class Save {
    public static String fileName = "system.ptcls";
    public static String controlSeparator = "-- Control Separator --";
    public static String particleSeparator = "-- Below are Particles --";
    public static String springSeparator = "-- Below are Springs --";
    static Particle[] particles;
    public static void save(ParticleSystem system) {
        try (PrintStream out  = new PrintStream(fileName)){
            out.println(controlSeparator);
            out.println(system.useGravity.getValue());
            out.println(system.gravity.getValue());
            out.println(system.viscousDamping.getValue());
            out.println(system.springStiffness.getValue());
            out.println(system.springDamping.getValue());
            out.println(system.iterations.getValue());
            out.println(system.restitution.getValue());

            out.println(particleSeparator);
            out.println(system.getParticles().size());
            for (Particle pt : system.getParticles()) {
                out.print(pt.index);
                out.print(";" + pt.p0.x);
                out.print(";" + pt.p0.y);
                out.print(";" + pt.v0.x);
                out.print(";" + pt.v0.y);
                out.print(";" + pt.p.x);
                out.print(";" + pt.p.y);
                out.print(";" + pt.v.x);
                out.print(";" + pt.v.y);
                out.print(";" + pt.f.x);
                out.print(";" + pt.f.y);
                out.print(";" + pt.pinned + "\n");
            }

            out.println(springSeparator);
            out.println(system.getSprings().size());
            for (Spring s : system.getSprings()) {
                out.print(s.p1.index);
                out.print(";" + s.p2.index);
                out.print(";" + s.l0 + "\n");
            }
        } catch(Exception ex){
            System.out.println("Save did not work");
        }
    }
    
    public static Particle createParticle(String[] particle_string) {
        Particle pt = new Particle(Double.valueOf(particle_string[1]), Double.valueOf(particle_string[2]), Double.valueOf(particle_string[3]), Double.valueOf(particle_string[4]) );
        pt.index = Integer.valueOf(particle_string[0]);
        pt.p.set(Double.valueOf(particle_string[5]), Double.valueOf(particle_string[6]));
        pt.v.set(Double.valueOf(particle_string[7]), Double.valueOf(particle_string[8]));
        pt.f.set(Double.valueOf(particle_string[9]), Double.valueOf(particle_string[10]));
        pt.pinned = Boolean.valueOf(particle_string[11]);
        return pt;
    }

    @SuppressWarnings("resource")
	public static void load(ParticleSystem system) {
        try {
            FileReader file_reader = new FileReader(fileName);
            BufferedReader reader = new BufferedReader(file_reader);
            system.clearParticles();
            if (reader.readLine().equals(controlSeparator)) {
                system.useGravity.setValue(Boolean.valueOf(reader.readLine()));
                system.gravity.setValue(Double.valueOf(reader.readLine()));
                system.viscousDamping.setValue(Double.valueOf(reader.readLine()));
                system.springStiffness.setValue(Double.valueOf(reader.readLine()));
                system.springDamping.setValue(Double.valueOf(reader.readLine()));
                system.iterations.setValue(Integer.valueOf(reader.readLine()));
                system.restitution.setValue(Double.valueOf(reader.readLine()));
            } else {
                throw new Exception("File Format Unknown");
            }
            if (reader.readLine().equals(particleSeparator)) {
                int n = Integer.valueOf(reader.readLine());
                particles = new Particle[n];
                for (int i = 0; i < n; i++) {
                    String[] particle_string = reader.readLine().split(";+");
                    Particle newpt = createParticle(particle_string);
                    system.getParticles().add(newpt);
                    particles[i] = newpt;
                }
            } else {
                throw new Exception("File Format Unknown");
            }
            if (reader.readLine().equals(springSeparator)) {
                int n = Integer.valueOf(reader.readLine());
                for (int i = 0; i < n; i++) {
                    String[] spring_string = reader.readLine().split(";+");
                    Particle tmp1 = particles[Integer.valueOf(spring_string[0])];
                    Particle tmp2 = particles[Integer.valueOf(spring_string[1])];
                    Spring s = new Spring(tmp1, tmp2);
                    s.l0 = Double.valueOf(spring_string[2]);
                    system.getSprings().add(s);
                }
            } else {
                throw new Exception("Unknown File Format");
            }
            reader.close();
            file_reader.close();
        } catch(Exception ex){
            System.out.println("Load Failed. Do you have \"system.ptcls\" in the working directory?");
        }
    }



}
