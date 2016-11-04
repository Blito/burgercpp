#ifndef MATERIAL_H
#define MATERIAL_H

enum class material_id
{
    GEL, AIR, FAT, BONE, BLOOD, VESSEL, LIVER, KIDNEY, SUPRARRENAL, GALLBLADDER, SKIN
};

struct material
{
    float impedance, attenuation, mu0, mu1, sigma;
};

struct EnumClassHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

#endif // MATERIAL_H
