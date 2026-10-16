#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>

namespace {
constexpr double A = 6378137.0;
constexpr double F = 1.0 / 298.257223563;
constexpr double E2 = F * (2.0 - F);
constexpr double PI = 3.14159265358979323846;
constexpr size_t MAX_CELLS = (128ULL * 1024 * 1024) / (3 * sizeof(double));
constexpr size_t MAX_STEPS = 1000000;
constexpr size_t MAX_TOTAL_STEPS = 10000000;

void require(bool ok, const std::string &message) {
    if (!ok) throw std::runtime_error(message);
}

void nc_check(int code, const std::string &what) {
    if (code != NC_NOERR) throw std::runtime_error(what + ": " + nc_strerror(code));
}

struct NcFile {
    int id = -1;
    explicit NcFile(const std::string &path) {
        nc_check(nc_open(path.c_str(), NC_NOWRITE, &id), "open weather cube");
    }
    ~NcFile() { if (id >= 0) nc_close(id); }
    NcFile(const NcFile &) = delete;
    NcFile &operator=(const NcFile &) = delete;
};

std::string text_attribute(int nc, int var, const char *name) {
    nc_type type;
    size_t len;
    nc_check(nc_inq_att(nc, var, name, &type, &len), std::string("required attribute ") + name);
    if (type == NC_STRING && len == 1) {
        char *text = nullptr;
        nc_check(nc_get_att_string(nc, var, name, &text), name);
        size_t length = text ? std::strlen(text) : 0;
        std::string value;
        if (text && length <= 256) value.assign(text, length);
        nc_free_string(1, &text);
        require(length <= 256, std::string("expected short text attribute: ") + name);
        return value;
    }
    require(type == NC_CHAR && len <= 256, std::string("expected short text attribute: ") + name);
    std::string value(len, '\0');
    nc_check(nc_get_att_text(nc, var, name, value.data()), name);
    return value;
}

void validate_epoch(const std::string &s) {
    require(s.size() == 20, "epoch_utc must be YYYY-MM-DDTHH:MM:SSZ");
    for (size_t i = 0; i < s.size(); ++i) {
        char separator = i == 4 || i == 7 ? '-' : i == 10 ? 'T' : i == 13 || i == 16 ? ':' : i == 19 ? 'Z' : '\0';
        require(separator ? s[i] == separator : s[i] >= '0' && s[i] <= '9', "invalid epoch_utc format");
    }
    int year = std::stoi(s.substr(0, 4)), month = std::stoi(s.substr(5, 2)), day = std::stoi(s.substr(8, 2));
    const int days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    require(year >= 1 && month >= 1 && month <= 12, "invalid epoch_utc calendar date");
    bool leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    require(day >= 1 && day <= days[month - 1] + (month == 2 && leap) &&
            std::stoi(s.substr(11, 2)) <= 23 && std::stoi(s.substr(14, 2)) <= 59 &&
            std::stoi(s.substr(17, 2)) <= 59, "invalid epoch_utc calendar date/time (leap seconds unsupported)");
}

bool numeric_type(nc_type t) {
    return t == NC_BYTE || t == NC_UBYTE || t == NC_SHORT || t == NC_USHORT ||
           t == NC_INT || t == NC_UINT || t == NC_INT64 || t == NC_UINT64 || t == NC_FLOAT || t == NC_DOUBLE;
}

double default_fill(nc_type t) {
    switch (t) {
        case NC_BYTE: return NC_FILL_BYTE;
        case NC_UBYTE: return NC_FILL_UBYTE;
        case NC_SHORT: return NC_FILL_SHORT;
        case NC_USHORT: return NC_FILL_USHORT;
        case NC_INT: return NC_FILL_INT;
        case NC_UINT: return NC_FILL_UINT;
        case NC_INT64: return static_cast<double>(NC_FILL_INT64);
        case NC_UINT64: return static_cast<double>(NC_FILL_UINT64);
        case NC_FLOAT: return NC_FILL_FLOAT;
        default: return NC_FILL_DOUBLE;
    }
}

std::vector<double> read_variable(int nc, const char *name, const std::vector<int> &dims,
                                  size_t count, const char *units) {
    int var, ndims;
    nc_type type;
    nc_check(nc_inq_varid(nc, name, &var), name);
    nc_check(nc_inq_var(nc, var, nullptr, &type, &ndims, nullptr, nullptr), name);
    require(numeric_type(type), std::string(name) + " must be numeric");
    require(ndims == static_cast<int>(dims.size()), std::string(name) + " has wrong dimension count");
    std::vector<int> actual(dims.size());
    nc_check(nc_inq_vardimid(nc, var, actual.data()), name);
    require(actual == dims, std::string(name) + " has wrong dimensions/order");
    require(text_attribute(nc, var, "units") == units, std::string(name) + " has wrong units; expected " + units);
    for (const char *attr : {"scale_factor", "add_offset"}) {
        int status = nc_inq_att(nc, var, attr, nullptr, nullptr);
        require(status == NC_ENOTATT, std::string(name) + " packed input rejected: " + attr);
    }
    std::vector<double> fills;
    for (const char *attr : {"_FillValue", "missing_value"}) {
        nc_type atype;
        size_t len;
        bool is_fill = std::string(attr) == "_FillValue";
        int status = nc_inq_att(nc, var, attr, &atype, &len);
        if (status == NC_ENOTATT) {
            if (is_fill) fills.push_back(default_fill(type));
            continue;
        }
        nc_check(status, attr);
        require(numeric_type(atype) && len >= 1 && len <= (is_fill ? 1U : 64U), std::string("invalid ") + attr);
        size_t old = fills.size();
        fills.resize(old + len);
        nc_check(nc_get_att_double(nc, var, attr, fills.data() + old), attr);
    }
    std::vector<double> values(count);
    nc_check(nc_get_var_double(nc, var, values.data()), name);
    for (double v : values) {
        require(std::isfinite(v) && std::find(fills.begin(), fills.end(), v) == fills.end(),
                std::string(name) + " contains nonfinite or missing/fill values");
    }
    return values;
}

struct Cell { size_t i; double t; };

Cell locate(const std::vector<double> &axis, double x, double tolerance, const char *name) {
    require(std::isfinite(x) && x >= axis.front() - tolerance && x <= axis.back() + tolerance,
            std::string("ray outside model ") + name + " coverage; no extrapolation");
    x = std::clamp(x, axis.front(), axis.back());
    auto it = std::upper_bound(axis.begin(), axis.end(), x);
    size_t i = it == axis.end() ? axis.size() - 2 : static_cast<size_t>(it - axis.begin() - 1);
    return {i, (x - axis[i]) / (axis[i + 1] - axis[i])};
}

struct Cube {
    std::vector<double> height, latitude, longitude, pressure, temperature, humidity;
    std::string epoch;

    explicit Cube(const std::string &path) {
        NcFile file(path);
        require(text_attribute(file.id, NC_GLOBAL, "height_reference") == "WGS84_ellipsoid",
                "height_reference must be WGS84_ellipsoid");
        epoch = text_attribute(file.id, NC_GLOBAL, "epoch_utc");
        validate_epoch(epoch);
        const char *names[] = {"height", "latitude", "longitude"};
        std::vector<int> dims(3);
        size_t lengths[3], cells = 1;
        for (size_t i = 0; i < 3; ++i) {
            nc_check(nc_inq_dimid(file.id, names[i], &dims[i]), names[i]);
            nc_check(nc_inq_dimlen(file.id, dims[i], &lengths[i]), names[i]);
            require(lengths[i] >= 2, std::string(names[i]) + " needs at least two nodes");
            require(lengths[i] <= MAX_CELLS / cells, "cube exceeds 128 MiB decoded three-variable limit");
            cells *= lengths[i];
        }
        height = read_variable(file.id, "height", {dims[0]}, lengths[0], "m");
        latitude = read_variable(file.id, "latitude", {dims[1]}, lengths[1], "degrees_north");
        longitude = read_variable(file.id, "longitude", {dims[2]}, lengths[2], "degrees_east");
        for (const auto *axis : {&height, &latitude, &longitude})
            for (size_t i = 1; i < axis->size(); ++i)
                require((*axis)[i] > (*axis)[i - 1], "coordinates must be strictly increasing; no flipping/wrapping");
        require(latitude.front() >= -90 && latitude.back() <= 90, "latitude must be within [-90,90]");
        require(longitude.front() >= -180 && longitude.back() <= 180, "longitude must be within [-180,180]; no wrapping");
        require(height.front() > -A * (1 - E2), "height lies below valid WGS84 normal coordinates");
        pressure = read_variable(file.id, "pressure", dims, cells, "Pa");
        temperature = read_variable(file.id, "temperature", dims, cells, "K");
        humidity = read_variable(file.id, "specific_humidity", dims, cells, "kg kg-1");
        for (size_t i = 0; i < cells; ++i)
            require(pressure[i] > 0 && temperature[i] > 0 && humidity[i] >= 0 && humidity[i] < 1,
                    "cube requires p>0, T>0 and 0<=q<1 at every node");
    }

    double refractivity(double lon, double lat, double h) const {
        // Tolerances only absorb ECEF conversion roundoff at boundary nodes.
        Cell z = locate(height, h, 1e-7, "height");
        Cell y = locate(latitude, lat, 1e-10, "latitude");
        Cell x = locate(longitude, lon, 1e-10, "longitude");
        auto interp = [&](const std::vector<double> &v) {
            double sum = 0;
            for (size_t k = 0; k < 2; ++k)
                for (size_t j = 0; j < 2; ++j)
                    for (size_t i = 0; i < 2; ++i) {
                        size_t index = ((z.i + k) * latitude.size() + y.i + j) * longitude.size() + x.i + i;
                        sum += v[index] * (k ? z.t : 1 - z.t) * (j ? y.t : 1 - y.t) * (i ? x.t : 1 - x.t);
                    }
            return sum;
        };
        double p = interp(pressure) * 0.01, t = interp(temperature), q = interp(humidity);
        double e = q * p / (0.622 + (1 - 0.622) * q);
        double n = 77.6 * ((p - e) / t) + 70.4 * (e / t) + 3.739e5 * (e / t) / t;
        require(std::isfinite(n) && n > 0, "nonfinite/nonpositive refractivity");
        return n;
    }
};

using Vec = std::array<double, 3>;

Vec to_ecef(double lon, double lat, double h) {
    lon *= PI / 180; lat *= PI / 180;
    double n = A / std::sqrt(1 - E2 * std::sin(lat) * std::sin(lat));
    return {(n + h) * std::cos(lat) * std::cos(lon), (n + h) * std::cos(lat) * std::sin(lon),
            (n * (1 - E2) + h) * std::sin(lat)};
}

Vec to_geodetic(const Vec &p) {
    double r = std::hypot(p[0], p[1]);
    double lat = std::atan2(p[2], r * (1 - E2));
    bool converged = false;
    for (int i = 0; i < 64; ++i) {
        double s = std::sin(lat), n = A / std::sqrt(1 - E2 * s * s);
        double next = std::atan2(p[2] + E2 * n * s, r);
        if (std::abs(next - lat) <= 2e-15) { lat = next; converged = true; break; }
        lat = next;
    }
    require(converged, "ECEF to WGS84 conversion failed to converge");
    double s = std::sin(lat);
    double h = r * std::cos(lat) + p[2] * s - A * std::sqrt(1 - E2 * s * s);
    require(std::isfinite(h), "nonfinite ECEF position/height");
    return {std::atan2(p[1], p[0]) * 180 / PI, lat * 180 / PI, h};
}

double number(const std::string &token) {
    char *end = nullptr;
    errno = 0;
    double v = std::strtod(token.c_str(), &end);
    require(end != token.c_str() && *end == '\0' && errno != ERANGE && std::isfinite(v),
            "expected a finite numeric value: " + token);
    return v;
}

struct Ray {
    std::array<double, 6> values;
    Vec origin, direction;
    explicit Ray(const std::array<double, 6> &v, const Cube &cube) : values(v) {
        require(v[0] >= -180 && v[0] <= 180 && v[1] >= -90 && v[1] <= 90, "invalid ray longitude/latitude");
        require(v[2] >= cube.height.front() && v[2] <= cube.height.back(), "ray start height outside model");
        require(v[0] >= cube.longitude.front() && v[0] <= cube.longitude.back() &&
                v[1] >= cube.latitude.front() && v[1] <= cube.latitude.back(), "ray start outside horizontal model coverage");
        double norm = std::hypot(v[3], v[4], v[5]);
        require(std::isfinite(norm) && std::abs(norm - 1) <= 1e-6 && v[5] > 0,
                "LOS ENU must be unit length (tolerance 1e-6) with up>0");
        origin = to_ecef(v[0], v[1], v[2]);
        double lon = v[0] * PI / 180, lat = v[1] * PI / 180;
        double e = v[3] / norm, n = v[4] / norm, u = v[5] / norm;
        direction = {-e * std::sin(lon) - n * std::sin(lat) * std::cos(lon) + u * std::cos(lat) * std::cos(lon),
                     e * std::cos(lon) - n * std::sin(lat) * std::sin(lon) + u * std::cos(lat) * std::sin(lon),
                     n * std::cos(lat) + u * std::sin(lat)};
    }

    Vec at(double s) const {
        return to_geodetic({origin[0] + s * direction[0], origin[1] + s * direction[1], origin[2] + s * direction[2]});
    }

    double length(double top) const {
        if (values[2] == top) return 0;
        double lo = 0, hi = std::max(1.0, top - values[2]);
        int expansions = 0;
        while (at(hi)[2] < top) {
            require(++expansions <= 128 && hi < std::numeric_limits<double>::max() / 2,
                    "cannot bracket top-height intersection");
            hi *= 2;
        }
        for (int i = 0; i < 100; ++i) {
            double mid = lo + (hi - lo) / 2;
            if (mid == lo || mid == hi) break;
            if (at(mid)[2] <= top) lo = mid; else hi = mid;
        }
        return lo;
    }
};

std::array<double, 2> integrate(const Cube &cube, const Ray &ray, double step, size_t &total) {
    double length = ray.length(cube.height.back());
    if (length == 0) return {0, 0};
    double needed = std::ceil(length / step);
    require(std::isfinite(needed) && needed <= MAX_STEPS, "ray exceeds 1000000 integration steps; increase step_m");
    size_t steps = std::max<size_t>(1, static_cast<size_t>(needed));
    require(steps <= MAX_TOTAL_STEPS - total, "run exceeds 10000000 total integration steps");
    total += steps;
    double ds = length / steps;
    double previous = cube.refractivity(ray.values[0], ray.values[1], ray.values[2]);
    double sum = 0, compensation = 0;
    for (size_t i = 1; i <= steps; ++i) {
        Vec point = ray.at(i == steps ? length : ds * i);
        require(std::abs(point[0] - ray.values[0]) <= 180, "ray crosses longitude seam; no wrapping");
        double current = cube.refractivity(point[0], point[1], point[2]);
        double term = (0.5 * previous + 0.5 * current) * 1e-6 * ds - compensation;
        double next = sum + term;
        compensation = (next - sum) - term;
        sum = next;
        previous = current;
    }
    require(std::isfinite(sum), "nonfinite integrated delay");
    return {sum, length};
}

struct Output {
    std::string temporary, destination;
    FILE *file = nullptr;
    explicit Output(const std::string &path) : destination(path) {
        struct stat info;
        int status = lstat(path.c_str(), &info);
        require(status != 0 && errno == ENOENT, "output already exists or cannot be inspected (inputs must differ from output)");
        std::string pattern = path + ".tmp.XXXXXX";
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        int fd = mkstemp(buffer.data());
        require(fd >= 0, "cannot create temporary output: " + std::string(std::strerror(errno)));
        temporary = buffer.data();
        file = fdopen(fd, "w");
        if (!file) {
            close(fd); unlink(temporary.c_str());
            throw std::runtime_error("cannot open temporary output stream");
        }
    }
    ~Output() {
        if (file) fclose(file);
        if (!temporary.empty()) unlink(temporary.c_str());
    }
    Output(const Output &) = delete;
    Output &operator=(const Output &) = delete;
    void publish() {
        require(fflush(file) == 0 && !ferror(file), "write output failed");
        require(fsync(fileno(file)) == 0, "sync output failed");
        int status = fclose(file);
        file = nullptr;
        require(status == 0, "close output failed");
        // Hard-link publication is atomic and never replaces an existing destination.
        require(link(temporary.c_str(), destination.c_str()) == 0,
                "publish output failed (will not overwrite): " + std::string(std::strerror(errno)));
    }
};

void help() {
    std::cout <<
        "Usage: atm_raytrace weather.nc rays.txt delays.txt [step_m]\n"
        "Independent straight-ray LOS excess path length, NOT a raw WRF/ERA5 reader.\n"
        "Input must be a normalized netCDF cube at ONE UTC epoch, not raw WRF/ERA5:\n"
        "  height(height) [m, WGS84 ellipsoid], latitude(latitude) [degrees_north],\n"
        "  longitude(longitude) [degrees_east]: finite, strictly increasing 1D, >=2 nodes.\n"
        "  latitude in [-90,90], longitude in [-180,180]; no flipping or seam wrapping.\n"
        "  Heights must exceed -a*(1-e^2), the WGS84 normal-coordinate lower bound.\n"
        "  pressure(height,latitude,longitude) [Pa], temperature(same dimensions) [K],\n"
        "  specific_humidity(same dimensions) [kg kg-1]; exact names/order/units.\n"
        "  All nodes: finite p>0, T>0, 0<=q<1; no missing/fill values.\n"
        "  scale_factor/add_offset attributes are rejected: decode before use.\n"
        "  Global text attributes: height_reference='WGS84_ellipsoid',\n"
        "  epoch_utc='YYYY-MM-DDTHH:MM:SSZ' (valid calendar, no leap seconds).\n"
        "  Three fields decoded as doubles must total <=128 MiB. No time interpolation.\n"
        "Rays: six finite columns lon_deg lat_deg ellipsoid_height_m east north up.\n"
        "  # comments allowed; ENU norm within 1e-6 of 1 (normalized), up>0.\n"
        "  Start and every integration node must be in cube coverage; no extrapolation.\n"
        "  Boundary roundoff tolerance only: 1e-7 m height, 1e-10 degrees horizontal.\n"
        "step_m: positive finite maximum segment length, default 100 m; at most\n"
        "  1000000 steps/ray and 10000000 steps/run. WGS84 ECEF geometry retains curvature.\n"
        "Output: lon lat h E N U truncated_los_delay_m path_length_m, in ray order.\n"
        "  Trapezoidal integral of refractivity along the straight LOS from start\n"
        "  to model top height only; NO ray bending or delay above model top.\n"
        "  Not phase/displacement; no automatic main-pipeline integration or sign convention.\n"
        "Output must differ from inputs and must not exist; atomic no-overwrite publication.\n";
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
            help(); return 0;
        }
        if (argc != 4 && argc != 5) { help(); return 2; }
        double step = argc == 5 ? number(argv[4]) : 100.0;
        require(step > 0, "step_m must be positive and finite");
        Cube cube(argv[1]);
        std::ifstream input(argv[2]);
        require(input.is_open(), "cannot open rays input");
        input.imbue(std::locale::classic());
        Output output(argv[3]);
        std::fprintf(output.file, "# epoch_utc=%s step_m=%.17g height_reference=WGS84_ellipsoid\n", cube.epoch.c_str(), step);
        std::fprintf(output.file, "# truncated straight-ray LOS excess path length: start to model_top_height_m=%.17g only; no bending or above-top delay; not phase/displacement\n", cube.height.back());
        std::fprintf(output.file, "# lon lat h E N U truncated_los_delay_m path_length_m\n");
        std::string line;
        size_t line_number = 0, rays = 0, total_steps = 0;
        while (std::getline(input, line)) {
            ++line_number;
            try {
                require(line.size() <= 8192, "ray line too long");
                size_t comment = line.find('#');
                if (comment != std::string::npos) line.resize(comment);
                std::istringstream fields(line);
                fields.imbue(std::locale::classic());
                std::string token;
                if (!(fields >> token)) continue;
                std::array<double, 6> values;
                values[0] = number(token);
                for (size_t i = 1; i < values.size(); ++i) {
                    require(static_cast<bool>(fields >> token), "ray must have exactly six columns");
                    values[i] = number(token);
                }
                require(!(fields >> token), "ray must have exactly six columns");
                Ray ray(values, cube);
                auto result = integrate(cube, ray, step, total_steps);
                for (double v : values) std::fprintf(output.file, "%.17g ", v);
                std::fprintf(output.file, "%.17g %.17g\n", result[0], result[1]);
                require(!ferror(output.file), "write output failed");
                ++rays;
            } catch (const std::exception &e) {
                throw std::runtime_error("rays line " + std::to_string(line_number) + ": " + e.what());
            }
        }
        require(!input.bad(), "error reading rays input");
        require(rays > 0, "no rays in input");
        output.publish();
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "atm_raytrace: " << e.what() << '\n';
        return 1;
    }
}
